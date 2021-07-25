/********* Sequence input routines for CLUSTAL W *******************/
/* DES was here.  FEB. 1994 */
/* Now reads PILEUP/MSF and CLUSTAL alignment files */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "score.h"	

#define MIN(a,b) ((a)<(b)?(a):(b))



/*
*	Prototypes
*/

static void check_infile(sint *noseqs);
static sint count_msf_seqs(void);
static char * get_seq(char *sname,sint *len);
static char * get_msf_seq(char *sname,sint *len,sint seqno);


/*
 *	Global variables
 */
FILE *fin;
static sint debug;

static sint seqFormat;
static char *formatNames[] = {"unknown","EMBL/Swiss-Prot","PIR",
			      "Pearson","GDE","Clustal","Pileup/MSF","RSF","USER","PHYLIP"};

static char * get_msf_seq(char *sname,sint *len,sint seqno)
/* read the seqno_th. sequence from a PILEUP multiple alignment file */
{
	static char line[MAXLINE+1];
	char *seq = NULL;
	sint i,j,k;
	unsigned char c;

	fseek(fin,0,0); 		/* start at the beginning */

	*len=0;				/* initialise length to zero */
        for(i=0;;i++) {
		if(fgets(line,MAXLINE+1,fin)==NULL) return NULL; /* read the title*/
		if(linetype(line,"//") ) break;		    /* lines...ignore*/
	}

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!blankline(line)) {

			for(i=0;i<seqno;i++) fgets(line,MAXLINE+1,fin);
                        for(j=0;j<=strlen(line);j++) if(line[j] != ' ') break;
			for(k=j;k<=strlen(line);k++) if(line[k] == ' ') break;
			strncpy(sname,line+j,MIN(MAXNAMES,k-j)); 
			sname[MIN(MAXNAMES,k-j)]=EOS;
			rtrim(sname);
                       	blank_to_(sname);
			if(seq==NULL)
				seq=(char *)ckalloc((MAXLINE+2)*sizeof(char));
			else
				seq=(char *)ckrealloc(seq,((*len)+MAXLINE+2)*sizeof(char));
			for(i=k;i<MAXLINE;i++) {
				c=line[i];
				if(c == '.' || c == '~' ) c = GAP2;
				if(c == '*') c = 'X';
				if(c == '\n' || c == EOS) break; /* EOL */
				if(isalpha(c)||c==GAP2) seq[(*len)++]=c;
			}

			for(i=0;;i++) {
				if(fgets(line,MAXLINE+1,fin)==NULL) return seq;
				if(blankline(line)) break;
			}
		}
	}
	return seq;
}

static char * get_seq(char *sname,sint *len)
{
	static char line[MAXLINE+1];
	char *seq = NULL;
	sint i, j, l, ti, oi, g, offset = 0;
	sint i1,i2;
        unsigned char c=EOS;
	Boolean got_seq=FALSE;
	Boolean found;

	ti=0;
	switch(seqFormat) {

/************************************/
		case PEARSON:
			while(*line != '>')
				fgets(line,MAXLINE+1,fin);
			
                        for(i=1;i<=strlen(line);i++)  /* DES */
				if(line[i] != ' ') break;
			strncpy(sname,line+i,MAXNAMES); /* remember entryname */
                        for(i=1;i<=strlen(sname);i++)  /* DES */
				if(sname[i] == ' ') break;
			sname[i]=EOS;
			rtrim(sname);
                        blank_to_(sname);

			*len=0;
			while(fgets(line,MAXLINE+1,fin)) {
				if(seq==NULL)
					seq=(char *)ckalloc((MAXLINE+2)*sizeof(char));
				else
					seq=(char *)ckrealloc(seq,((*len)+MAXLINE+2)*sizeof(char));
				for(i=0;i<MAXLINE;i++) {
					c=line[i];
				if(c == '\n' || c == EOS || c == '>')
					break;			/* EOL */
			
				if(isalpha(c)||c==GAP2) seq[(*len)++]=c;
			}
			if(c == '>') break;
			}
		break;
/**********************************************/
	}
	
	seq[*len+1]=EOS;

	return seq;
}

sint readseqs(char *filename,ALNPTR mult_aln)
{
	static char *seq1,sname1[MAXNAMES+1];
	char tstr[MAXTITLES+1];
	char tmp[MAXTITLES+1];
	sint i,j,k;
	sint no_seqs;
	sint group;
	sint seq_type;
	static sint l1;
	static sint max_aln_length;
	static Boolean dnaflag1;
	Boolean found;
	
	if(*filename == EOS) return -1;
	if((fin=fopen(filename,"r"))==NULL) {
		error("Could not open sequence file %s",filename);
		return -1;      /* DES -1 => file not found */
	}
	no_seqs=0;

	check_infile(&no_seqs);

	if(no_seqs == 0)
		return 0;       /* return the number of seqs. (zero here)*/

	alloc_aln(no_seqs,mult_aln);
	for(i=0;i<no_seqs;i++) {    /* get the seqs now*/
		if(seqFormat == MSF){
		    	seq1=get_msf_seq(sname1,&l1,i);
		}
		else
			seq1=get_seq(sname1,&l1);

		if(seq1==NULL) break;
		if (l1 > max_aln_length) max_aln_length = l1;
		mult_aln->seqs[i].len=l1;                   /* store the length */
		strcpy(mult_aln->seqs[i].name,sname1);              /*    "   "  name   */

		alloc_seq(&mult_aln->seqs[i],l1);

		for(j=0;j<l1;j++) {
			if(isalpha(seq1[j])) mult_aln->seqs[i].data[j]=tolower(seq1[j]);
			else mult_aln->seqs[i].data[j]=GAP2;
		}
		if(seq1!=NULL) seq1=ckfree(seq1);

		mult_aln->seqs[i].reslen=0;
		for(j=0;j<mult_aln->seqs[i].len;j++)
			if(isalpha(mult_aln->seqs[i].data[j])) mult_aln->seqs[i].reslen++;
	}

	max_aln_length *= 2;
	for(i=0;i<no_seqs;i++)
	{
		if(mult_aln->seqs[i].len>max_aln_length)
			max_aln_length=mult_aln->seqs[i].len;
	}
	
	fclose(fin);
			
	mult_aln->nseqs=no_seqs;
	return no_seqs;    /* return the number of seqs. read in this call */
}

static void check_infile(sint *nseqs)
{
	char line[MAXLINE+1];
	sint i;	

	*nseqs=0;
	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!blankline(line)) 
			break;
	}

	for(i=strlen(line)-1;i>=0;i--)
		if(isgraph(line[i])) break;
	line[i+1]=EOS;
        
	for(i=0;i<=6;i++) line[i] = toupper(line[i]);

 	if( linetype(line,"PILEUP") ) {
		seqFormat = MSF;
	}
 	else if( strstr(line,"MSF") && line[strlen(line)-1]=='.' &&
                                 line[strlen(line)-2]=='.' ) {
		seqFormat = MSF;
	}
	else if(*line == '>') {	
		seqFormat=PEARSON; 
		(*nseqs)++;
	}
	else {
		seqFormat=UNKNOWN;
		return;
	}

	while(fgets(line,MAXLINE+1,fin) != NULL) {
		switch(seqFormat) {
			case PEARSON:
                                if( *line == '>' )
                                        (*nseqs)++;
                                break;
			case MSF:
				*nseqs = count_msf_seqs();
				fseek(fin,0,0);
				return;
			case USER:
			default:
				break;
		}
	}
	fseek(fin,0,0);
}


static sint count_msf_seqs(void)
{
/* count the number of sequences in a PILEUP alignment file */

	char line[MAXLINE+1];
	sint  nseqs;

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(linetype(line,"//")) break;
	}

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!blankline(line)) break;		/* Look for next non- */
	}						/* blank line */
	nseqs = 1;

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(blankline(line)) return nseqs;
		nseqs++;
	}

	return (sint)0;	/* if you got to here-funny format/no seqs.*/
}



