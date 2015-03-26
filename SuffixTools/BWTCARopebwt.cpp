//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTCARopebwt - Construct the BWT for a set of reads
// using Heng Li's ropebwt implementation
//
#include "BWTCARopebwt.h"
#include "bcr.h"
#include "SeqReader.h"
#include "StdAlnTools.h"
#include "BWTWriterBinary.h"
#include "BWTWriterAscii.h"
#include "SAWriter.h"

/*** ropebwt2 headers ROPEBWT2_VERSION r187 ***/
#include <zlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "rld0.h"
#include "rle.h"
#include "mrope.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define FLAG_FOR 0x1
#define FLAG_REV 0x2
#define FLAG_ODD 0x4
#define FLAG_BIN 0x8
#define FLAG_TREE 0x10
#define FLAG_THR 0x40
#define FLAG_LINE 0x100
#define FLAG_RLD 0x200
#define FLAG_NON 0x400
#define FLAG_CRLF 0x800
#define FLAG_CUTN 0x1000


static inline int kputsn(const char *p, int l, kstring_t *s)
{
	if (s->l + l + 1 >= s->m) {
		char *tmp;
		s->m = s->l + l + 2;
		kroundup32(s->m);
		if ((tmp = (char*)realloc(s->s, s->m))) s->s = tmp;
		else return EOF;
	}
	memcpy(s->s + s->l, p, l);
	s->l += l;
	s->s[s->l] = 0;
	return l;
}

static void liftrlimit() // increase the soft limit to hard limit
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	if (r.rlim_cur < r.rlim_max) r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}
double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}


static unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5
};

void BWTCA::runRopebwt(const std::string& input_filename, const std::string& bwt_out_name,
                       bool use_threads, bool do_reverse)
{
    // Initialize ropebwt
    std::string tmp_name = bwt_out_name + ".tmp";
    bcr_t* bcr = bcr_init(use_threads, tmp_name.c_str());

    size_t num_sequences = 0;
    size_t num_bases = 0;
    SeqReader reader(input_filename);
    SeqRecord record;
    while(reader.get(record))
    {
        if(do_reverse)
            record.seq.reverse();

        size_t l = record.seq.length();

        // Convert the string into the alphabet encoding expected by ropebwt
        uint8_t* s = new uint8_t[l];
        for(size_t i = 0; i < l; ++i) {
            char c = record.seq.get(i);
            s[i] = seq_nt6_table[(int)c];
        }

        // Send the sequence to ropebwt
        bcr_append(bcr, l, s);

        num_sequences += 1;
        num_bases += l;
        delete [] s;
    }

    // Build the BWT
    bcr_build(bcr);

    // write the BWT and SAI
    bcritr_t* itr = bcr_itr_init(bcr);
    const uint8_t* s;
    int l;

    BWTWriterBinary* out_bwt = new BWTWriterBinary(bwt_out_name);
    size_t num_symbols = num_bases + num_sequences;
    out_bwt->writeHeader(num_sequences, num_symbols, BWF_NOFMI);

    // Write each run
    while( (s = bcr_itr_next(itr, &l)) != 0 ) {
        for (int i = 0; i < l; ++i) {
            char c = "$ACGTN"[s[i]&7];
            int rl = s[i]>>3;
            for(int j = 0; j < rl; ++j)
                out_bwt->writeBWChar(c);
        }
    }

    free(itr);
    out_bwt->finalize();
    delete out_bwt;

    // Cleanup
    bcr_destroy(bcr);
}


void BWTCA::runRopebwt2(const std::string& input_filename, const std::string& bwt_out_name,
                       int thr_min, bool do_reverse)
{
	mrope_t *mr = 0;
	gzFile fp;
	kseq_t *ks;
	int64_t m = (int64_t)(.97 * 10 * 1024 * 1024 * 1024) + 1;;
	int i, block_len = ROPE_DEF_BLOCK_LEN, max_nodes = ROPE_DEF_MAX_NODES, so = MR_SO_IO;
	int flag = FLAG_FOR | FLAG_REV | FLAG_THR;
	kstring_t buf = { 0, 0, 0 };
	double ct, rt;

	liftrlimit();
	if (mr == 0) mr = mr_init(max_nodes, block_len, so);
	if (thr_min > 0) mr_thr_min(mr, thr_min);
	fp = gzopen( input_filename.c_str(), "rb");
	ks = kseq_init(fp);
	ct = cputime(); rt = realtime();

	for (;;) {
		int l;
		uint8_t *s;
		if (kseq_read(ks) < 0) break; // read fasta/fastq
		
		l = ks->seq.l;
		s = (uint8_t*)ks->seq.s;

		// change encoding according to seq_nt6_table
		for (i = 0; i < l; ++i) 
			s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;
				
		if(!do_reverse)
			for (i = 0; i < l>>1; ++i) { // reverse
				int tmp = s[l-1-i];
				s[l-1-i] = s[i]; s[i] = tmp;
			}
		
		//push into buffer, always process the forward strand
		if (m) kputsn((char*)ks->seq.s, ks->seq.l + 1, &buf);
		else mr_insert1(mr, s);

		//sort if buffer is full
		if (m && (int64_t) buf.l >= m) {
			mr_insert_multi(mr, buf.l, (uint8_t*)buf.s, flag&FLAG_THR);
			buf.l = 0;
		}		
	}

	if (m && buf.l) 
		mr_insert_multi(mr, buf.l, (uint8_t*)buf.s, flag&FLAG_THR);
	
	//Done BWT construction
	int64_t c[6];
	fprintf(stderr, "[%s] constructed FM-index in %.3f sec, %.3f CPU sec\n", __func__, realtime() - rt, cputime() - ct);
	mr_get_c(mr, c);
	fprintf(stderr, "[%s] symbol counts: ($, A, C, G, T, N) = (%ld, %ld, %ld, %ld, %ld, %ld)\n", __func__,
			(long)c[0], (long)c[1], (long)c[2], (long)c[3], (long)c[4], (long)c[5]);
			
	int64_t num_sequences = (long)c[0];
    int64_t num_symbols = (long)c[0]+(long)c[1]+(long)c[2]+(long)c[3]+(long)c[4]+(long)c[5];

	free(buf.s);
	kseq_destroy(ks);
	gzclose(fp);
	
	/*** output BWT char to BWTWriter ***/
    BWTWriterBinary* out_bwt = new BWTWriterBinary(bwt_out_name);
    out_bwt->writeHeader(num_sequences, num_symbols, BWF_NOFMI);

	mritr_t itr;
	const uint8_t *block;

	mr_itr_first(mr, &itr, 1);
	while ((block = mr_itr_next_block(&itr)) != 0) {
		const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
		while (q < end) {
			int c = 0;
			int64_t j, l;
			rle_dec1(q, c, l);
			for (j = 0; j < l; ++j) 
				//putchar("$ACGTN"[c]);
				out_bwt->writeBWChar("$ACGTN"[c]);
		}		
	}
	
    out_bwt->finalize();
    delete out_bwt;
}

