/*  bali_score_rv100.c  - version 1.00 					*/
/*                   							*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include "score.h"

typedef struct {
	int start_col;
	int end_col;
	int len;
	int nseqs;
	int *seq;
} SBLOCK,*SBLOCKPTR;

int checkref(FILE *fin);
void aln_seq(int i);
ALN read_msf(FILE *fin,int nseqs);
int countmsf(FILE *fin);
void code_refseq(ALNPTR ref_aln,ALNPTR test_aln,int refmaxlen,int **refseq_code,int *seq_xref,int maxlen,int **seq_code);
float md_score(ALNPTR mult_aln,SBLOCK block,float matrix[NUMRES][NUMRES]);
void get_default_matrix(void);


float matrix[NUMRES][NUMRES];

int main(int argc, char **argv)
{
	FILE *ifd,*tfd,*afd;
	int  err,i,j,k,l;
	int ix,n;
	int nseqs,refnseqs;
	int refmaxlen,maxlen,ix1,ix2;
	int format;
	int nblocks,len;
	int seqi,seqj;
	int start,end,start_col,end_col;
	int start1,start2;
	int test_start_col1,test_start_col2;
	int p1,p2,c1,c2;
	float max_sp_score,sp_score,max_col_score,col_score;
	float bmax_sp_score,bsp_score,block_sp_score,bmax_col_score,bcol_score,block_col_score;
        int **seq_code,**refseq_code;
	int tmp_col_score[1000];
	int tmp_maxcol_score[1000];
	char seq[MAXLINE+1];
	char method;
	Boolean eof,found;
	int *seq_xref;
	int pcid,pcseqs;
	float md;
	int meanpcid,ndisorder;
	ALN ref_aln;
	ALN test_aln;
	SBLOCK *block;

	if(argc!=3) {
		fprintf(stderr,"Usage: %s ref_aln test_aln \n",argv[0]);
		fprintf(stderr,"                where ref_aln       reference alignment in xml/msf format \n");
		fprintf(stderr,"                      test_aln      test alignment in msf format \n");
			return 1;
	}

/* open the reference aln file */

        if((ifd=fopen(argv[1],"r"))==NULL) {
            fprintf(stderr,"Cannot open reference aln file [%s]\n",argv[1]);
            return 1;
        }

/* read the reference alignment into refaln */
	format=checkref(ifd);

	if(format==1) {
		refnseqs=count_xml_seqs(ifd);
	}
	else {
            fprintf(stderr,"Bad format in reference aln file [%s]\n",argv[1]);
            return 1;
	}

        fseek(ifd,0,0);
	if(refnseqs==0) {
		fprintf(stderr,"Error: no sequences in %s\n",argv[1]);
		return 1;
	}
	ref_aln=read_xml(ifd,0);
	ref_aln.nseqs=refnseqs;
	
	refmaxlen=0;
	for(i=0;i<ref_aln.nseqs;i++)
		if(refmaxlen<ref_aln.seqs[i].len) refmaxlen=ref_aln.seqs[i].len;

/* read the test alignment into testaln */
	nseqs = readseqs(argv[2],&test_aln);
	if(nseqs==0) {
		fprintf(stderr,"Error: no sequences in %s\n",argv[2]);
		return 1;
	}

	if(nseqs != refnseqs) {
		fprintf(stderr,"Error: %d sequences in %s and %d in %s",refnseqs,argv[1],nseqs,argv[2]);
		return 1;
	}

	maxlen=0;
	for(i=0;i<test_aln.nseqs;i++)
		if(maxlen<test_aln.seqs[i].len) maxlen=test_aln.seqs[i].len;



/* cross-reference the sequence refnames, in case they're not in the same order
in the reference and test alignment files */
	seq_xref=(int *)ckalloc((refnseqs+2)*sizeof(int));
	for(i=0;i<refnseqs;i++)
		seq_xref[i]=-1;
	for(i=0;i<refnseqs;i++) {
		found=FALSE;
		for(j=0;j<test_aln.nseqs;j++) {
			if(strcasecmp(test_aln.seqs[j].name,ref_aln.seqs[i].name)==0)
			{
				found=TRUE;
				seq_xref[i]=j;
				break;
			}
		}
		if(found==FALSE) {
			fprintf(stderr,"Error: sequence %s not found in test aln %s\n",ref_aln.seqs[i].name,argv[2]);
			return 1;
		}
	}

	get_default_matrix();

        fprintf(stdout,"\nComparing test alignment in %s\nwith reference alignment in %s\n",argv[2],argv[1]);

/* count the core blocks */
	nblocks=0;
	for(seqi=0;seqi<refnseqs;seqi++) {
		for(i=0;i<ref_aln.ft[seqi].nentries[COREBLOCK];i++) {
			if(strcmp(ref_aln.ft[seqi].data[COREBLOCK][i].name,"LBLOCK")==0) continue;
			if(ref_aln.ft[seqi].data[COREBLOCK][i].color>nblocks) nblocks=ref_aln.ft[seqi].data[COREBLOCK][i].color;
		}
	}
	nblocks++;

/* get the core blocks */
	block=(SBLOCK *)ckalloc((nblocks+1)*sizeof(SBLOCK));
	for (i=0;i<nblocks;i++)
		block[i].seq=(int *)ckalloc((refnseqs+1)*sizeof(int));

	for(seqi=0;seqi<refnseqs;seqi++) {
		if(ref_aln.ft[seqi].nentries[SEQERRBLOCK]>0) continue;
		for(i=0;i<ref_aln.ft[seqi].nentries[COREBLOCK];i++) {
			if(strcmp(ref_aln.ft[seqi].data[COREBLOCK][i].name,"LBLOCK")==0) continue;
			start=ref_aln.ft[seqi].data[COREBLOCK][i].start;
			end=ref_aln.ft[seqi].data[COREBLOCK][i].end;
			k=ref_aln.ft[seqi].data[COREBLOCK][i].color;
			pos2col(ref_aln.seqs[seqi].data,start,end,&start_col,&end_col);
			if(block[k].nseqs==0 || block[k].start_col<start_col) block[k].start_col=start_col;
			if(block[k].nseqs==0 || block[k].end_col>end_col) block[k].end_col=end_col;
			block[k].len=end-start+1;
			block[k].seq[block[k].nseqs]=seqi;
			block[k].nseqs++;
		}
	}

        for(i=0;i<ref_aln.nseqs;i++) {
                for(j=0;j<MAXFTTYPE;j++) {
                        if(j!=DISORDER) continue;
                        for(k=0;k<ref_aln.ft[i].nentries[j];k++) {
                                pos2col(ref_aln.seqs[i].data,ref_aln.ft[i].data[j][k].start,ref_aln.ft[i].data[j][k].end,&ref_aln.ft[i].data[j][k].start_col,&ref_aln.ft[i].data[j][k].end_col);
                        }
                }
        }

/* code the reference alignment - assign to each residue the number of the column it's in
   gap positions are coded -1 */

        refseq_code=(int **)ckalloc((ref_aln.nseqs+2)*sizeof(int *));
        seq_code=(int **)ckalloc((ref_aln.nseqs+2)*sizeof(int *));
        for(i=0;i<ref_aln.nseqs;i++) {
                refseq_code[i]=(int *)ckalloc((refmaxlen+2)*sizeof(int));
                seq_code[i]=(int *)ckalloc((refmaxlen+2)*sizeof(int));
	}
        code_refseq(&ref_aln,&test_aln,refmaxlen,refseq_code,seq_xref,maxlen,seq_code);


/* calculate the score based on core blocks only */

	max_sp_score=sp_score=max_col_score=col_score=0;
	for(n=0;n<nblocks;n++) {
		if(block[n].nseqs<=1) continue;
		if(block[n].end_col<=block[n].start_col) continue;
		pcseqs=(float)block[n].nseqs*100.0/(float)ref_aln.nseqs;
		bmax_sp_score=bsp_score=0;
		start_col=block[n].start_col;
		end_col=block[n].end_col;
		for(k=start_col;k<=end_col;k++) {
			tmp_maxcol_score[k-start_col]=0;
			tmp_col_score[k-start_col]=0;
		}
		meanpcid=ix=0;
		for(i=0;i<block[n].nseqs;i++) {
			seqi=block[n].seq[i];
			for(j=i+1;j<block[n].nseqs;j++) {
				seqj=block[n].seq[j];
/* calculate %id of block in two sequences */
				len=pcid=0;
				for(k=start_col;k<=end_col;k++) {
					if(refseq_code[seqi][k]>=0 && refseq_code[seqj][k]>=0) {
						len++;
						if(ref_aln.seqs[seqi].data[k] == ref_aln.seqs[seqj].data[k])
							pcid++;
					}
				}
				if(len>0) pcid=(float)pcid*100.0/(float)len;
				else pcid=0;
				meanpcid+=pcid;
				ix++;
				for(k=start_col;k<=end_col;k++) {
					if(refseq_code[seqi][k]>=0 && refseq_code[seqi][k]==refseq_code[seqj][k]) {
						bmax_sp_score++;
						tmp_maxcol_score[k-start_col]++;
					}
					if(seq_code[seqi][k]>=0 && seq_code[seqi][k]==seq_code[seqj][k]) {
						bsp_score++;
						tmp_col_score[k-start_col]++;
					}
				}
			}
		}
		if(ix>0) meanpcid=(float)meanpcid/(float)ix;
		if(bmax_sp_score>0) block_sp_score=bsp_score/bmax_sp_score;
		else block_sp_score=0;

/* score the columns: 0 for a column misaligned, no. of seqs in block for a column aligned */
		bcol_score=bmax_col_score=0;
		for(k=start_col;k<=end_col;k++) {
			if(tmp_maxcol_score[k-start_col]>0) {
				/*bmax_col_score++;
				if(tmp_col_score[k-start_col]==tmp_maxcol_score[k-start_col]) bcol_score++;*/
				bmax_col_score+=block[n].nseqs;
				if(tmp_col_score[k-start_col]==tmp_maxcol_score[k-start_col]) bcol_score+=block[n].nseqs;
			}
		}
		if(bmax_col_score>0) block_col_score=bcol_score/bmax_col_score;
		else block_col_score=0;

/* calculate md score for the block */
		md=md_score(&ref_aln,block[n],matrix);

/* calculate percentage of block prdicted to be disordered */

                ndisorder=0;
                j=DISORDER;
                for(i=0;i<block[n].nseqs;i++) {
                        seqi=block[n].seq[i];
                        for(k=0;k<ref_aln.ft[seqi].nentries[j];k++) {
                                if((l=overlap(block[n].start_col,block[n].end_col,ref_aln.ft[seqi].data[j][k].start_col,ref_aln.ft[seqi].data[j][k].end_col))>0) {
                                        ndisorder+=l;
                                        break;
                                }
                        }
                }
		ndisorder=(int)((float)ndisorder*100.0/(float)(block[n].nseqs*block[n].len));


		max_sp_score+=bmax_sp_score;
		sp_score+=bsp_score;
		max_col_score+=bmax_col_score;
		col_score+=bcol_score;
	}

	if(max_sp_score>0) sp_score=sp_score/max_sp_score;
	else sp_score=0;
	if(max_col_score>0) col_score=col_score/max_col_score;
	else col_score=0;

        fprintf(stdout,"\n\tCS score= %.3f\n",col_score); 

	exit(0);

	
}

float md_score(ALNPTR mult_aln,SBLOCK block,float matrix[NUMRES][NUMRES])
{
        char c,c1;
        int i,j,k,s,s1,p,r;
        int len,is,ie;
        int nseqs;
        float mean,n,score;
        float dist,diff,sum;
        float seqvector[26],seqvector1[26];
        float resdist[26][26];
        float resfreq[26][26];

        nseqs=block.nseqs;

        for(i=0;i<26;i++) {
                for (r=0;r<26; r++)
                        seqvector[r]=matrix[r][i];
                resdist[i][i]=0.0;
                for(j=i+1;j<26;j++) {
                        for (r=0;r<26; r++)
                                seqvector1[r]=matrix[r][j];
                        resdist[i][j]=0.0;
                        for(r=0;r<26;r++) {
                                diff=seqvector1[r]-seqvector[r];
                                resdist[i][j]+=diff*diff;
                        }
                        resdist[i][j]=sqrt((double)resdist[i][j]);
                        resdist[j][i]=resdist[i][j];
                }
        }
        sum=0.0;
        for(p=block.start_col; p<=block.end_col; p++)
        {
                for(i=0;i<26;i++)
                        for(j=0;j<26;j++)
                                resfreq[i][j]=0.0;
/* calculate mean of seq distances */
                mean=0.0;
                for(i=0; i<nseqs; i++) {
                        s=block.seq[i];
                        if(isalpha(mult_aln->seqs[s].data[p])) {
                                c=toupper(mult_aln->seqs[s].data[p]);
                                for(j=i+1; j<nseqs; j++) {
                                        s1=block.seq[j];
                                        if(isalpha(mult_aln->seqs[s1].data[p])) {
                                                c1=toupper(mult_aln->seqs[s1].data[p]);
                                                resfreq[c-'A'][c1-'A']++;
                                        }
                                }
                        }
                }
                for(i=0;i<26;i++) {
                        for(j=0;j<26;j++) {
                                mean+=resdist[i][j]*resfreq[i][j];
                        }
                }
                mean/=(float)(nseqs)*(float)(nseqs-1)/2.0;

/* normalise score between 0 and 1 (where 1 is an identical column) */

                score=exp((double)(-mean)/(double)2.0);

                sum+=score;
        }

        return sum/(float)(block.end_col-block.start_col+1);;
}


void code_refseq(ALNPTR ref_aln,ALNPTR test_aln,int refmaxlen,int **refseq_code,int *seq_xref,int maxlen,int **seq_code)
{
        int seq,i,j;

/* assign column no. in reference alignment to each residue */

        for(seq=0;seq<ref_aln->nseqs;seq++) {
		j=0; /* column in test alignment */
                for(i=0;i<refmaxlen;i++) {
                        if(i>=ref_aln->seqs[seq].len) {
                                refseq_code[seq][i]=-1;
                        }
                        else if(ref_aln->seqs[seq].data[i]=='-') {
                                refseq_code[seq][i]=-1;
                        }
                        else {
                                refseq_code[seq][i]=i;
/* get the equivalent residue in the test alignment and assign column no. in test alignment */
				for(;j<maxlen;j++) 
					if(test_aln->seqs[seq_xref[seq]].data[j]!='-') break;
                                seq_code[seq][i]=j;
				j++;
                        }
                }
        }
}

int checkref(FILE *fin)
{

        char line[MAXLINE+1];

	int format;

	format=0;

        while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(strncmp(line,"<?xml",5)==0) {
			fprintf(stdout,"Using XML format reference alignment\n");
			format=1;
			break;
		}
        }
        fseek(fin,0,0);

	return format;
}

void fatal( char *msg,...)
{
        va_list ap;

        va_start(ap,msg);
        fprintf(stdout,"\n\nFATAL ERROR: ");
        vfprintf(stdout,msg,ap);
        fprintf(stdout,"\n\n");
        va_end(ap);
        exit(1);
}

void error( char *msg,...)
{
        va_list ap;

        va_start(ap,msg);
        fprintf(stdout,"\n\nERROR: ");
        vfprintf(stdout,msg,ap);
        fprintf(stdout,"\n\n");
        va_end(ap);
}

void warning( char *msg,...)
{
        va_list ap;

        va_start(ap,msg);
        fprintf(stdout,"\n\nWARNING: ");
        vfprintf(stdout,msg,ap);
        fprintf(stdout,"\n\n");
        va_end(ap);
}

void info( char *msg,...)
{
        va_list ap;

        va_start(ap,msg);
        fprintf(stdout,"\n");
        vfprintf(stdout,msg,ap);
        va_end(ap);
}


short gon250mt[]={
  24,
   0,   0,
   5,   0, 115,
  -3,   0, -32,  47,
   0,   0, -30,  27,  36,
 -23,   0,  -8, -45, -39,  70,
   5,   0, -20,   1,  -8, -52,  66,
  -8,   0, -13,   4,   4,  -1, -14,  60, 
  -8,   0, -11, -38, -27,  10, -45, -22,  40,
  -4,   0, -28,   5,  12, -33, -11,   6, -21,  32,
 -12,   0, -15, -40, -28,  20, -44, -19,  28, -21,  40,
  -7,   0,  -9, -30, -20,  16, -35, -13,  25, -14,  28,  43,
  -3,   0, -18,  22,   9, -31,   4,  12, -28,   8, -30, -22,  38,
   3,   0, -31,  -7,  -5, -38, -16, -11, -26,  -6, -23, -24,  -9,  76,
  -2,   0, -24,   9,  17, -26, -10,  12, -19,  15, -16, -10,   7,  -2,  27,
  -6,   0, -22,  -3,   4, -32, -10,   6, -24,  27, -22, -17,   3,  -9,  15,  47,
  11,   0,   1,   5,   2, -28,   4,  -2, -18,   1, -21, -14,   9,   4,   2,  -2,  22,
   6,   0,  -5,   0,  -1, -22, -11,  -3,  -6,   1, -13,  -6,   5,   1,   0,  -2,  15,  25,
   1,   0,   0, -29, -19,   1, -33, -20,  31, -17,  18,  16, -22, -18, -15, -20, -10,   0,  34,
 -36,   0, -10, -52, -43,  36, -40,  -8, -18, -35,  -7, -10, -36, -50, -27, -16, -33, -35, -26, 142,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -22,   0,  -5, -28, -27,  51, -40,  22,  -7, -21,   0,  -2, -14, -31, -17, -18, -19, -19, -11,  41,   0,  78,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

void get_default_matrix(void)
{
   int gg_score = 0;
   int gr_score = 0;
   int i, j, k, ix = 0;
   int ti, tj;
   int  maxres;
   int av1,av2,av3,min, max;
   int max_aa;
   float scale=0.1;
        char *amino_acid_order = "ABCDEFGHIKLMNPQRSTVWXYZ";
   char *amino_acid_codes   =    "ABCDEFGHIJKLMNOPQRSTUVWXYZ-";
   static short  *xref, *matptr;
   short    def_aa_xref[NUMRES+1];

   char c1,c2;

   max_aa = strlen(amino_acid_codes)-2;

/*
   set up cross-reference for default matrices hard-coded in matrices.h
*/
   for (i=0;i<NUMRES;i++) def_aa_xref[i] = -1;

   maxres = 0;
   for (i=0;(c1=amino_acid_order[i]);i++)
     {
         for (j=0;(c2=amino_acid_codes[j]);j++)
          {
           if (c1 == c2)
               {
                  def_aa_xref[i] = j;
                  maxres++;
                  break;
               }
          }
     }



   matptr = gon250mt;
   xref=def_aa_xref;
/*
   default - set all scores to 0
*/
   for (i=0;i<=max_aa;i++)
      for (j=0;j<=max_aa;j++)
          matrix[i][j] = 0;

   ix = 0;
   maxres = 0;
   for (i=0;i<=max_aa;i++)
    {
      ti = xref[i];
      for (j=0;j<=i;j++)
       {
          tj = xref[j];
          if ((ti != -1) && (tj != -1))
            {
               k = matptr[ix];
               if (ti==tj)
                  {
                     matrix[ti][ti] = k * scale;
                     maxres++;
                  }
               else
                  {
                     matrix[ti][tj] = k * scale;
                     matrix[tj][ti] = k * scale;
                  }
               ix++;
            }
       }
    }

   --maxres;

   av1 = av2 = av3 = 0;
   for (i=0;i<=max_aa;i++)
    {
      for (j=0;j<=i;j++)
       {
           av1 += matrix[i][j];
           if (i==j)
              {
                 av2 += matrix[i][j];
              }
           else
              {
                 av3 += matrix[i][j];
              }
       }
    }

}

