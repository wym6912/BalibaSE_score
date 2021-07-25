#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "score.h"
#include <expat.h>

#define BUFFSIZE	8192 	/* size of buffer for reading xml file */

char Buff[BUFFSIZE]; 	/* buffer for reading xml file */

int Depth;		/* depth of xml tree (not currently used!) */
char *Element;		/* element name */
char *Content;		/* element contents */
char **Attributes;	/* list of attributes */
int  NAttributes;	/* number of attributes */
int Nseqs;		/* number of sequences */
float Alnscore;		/* alignment score */
int Nftable[MAXFTTYPE];	/* number of entries in feature table */
int Ftype;		/* feature type */
int Ngorefs;		/* number of GO cross-references */
int Ntaxons;		/* number of taxon in lineage */
int Nkeywords;		/* number of keywords */
int Ncolscores;		/* number of column scores */
int Fseq;		/* first sequence - start index to Maln */
int Cresidue1;		/* contact residue 1 */
int NCresidue;		/* number of contact residues */
int NC;
RCONTACT *Rcontacts;
ALN Maln;

#define START 0
#define END 1

static void do_element(void);
int getfloatargs(char *inline1,float *args,int max);

/* This is the application specific stuff. When we get here, the element name is contained in
Element, any attributes are in NAttributes and Attributes, the data is in Content */

static void do_element(void)
{
	int i,j,k,l,n;
	int cpos1,cpos2;
	char ec_text[30];
	int ec[4];

	if(strcmp(Element,"sequence")==0) {
/* this is the start of a new sequence, so reset the data for this sequence */
		for(i=0;i<MAXFTTYPE;i++)
			Nftable[i]=(-1);
		Ngorefs=(-1);
		Nkeywords=(-1);
		Ntaxons=(-1);
		Ncolscores=(-1);
		NC=0;
		if(strcasecmp(Attributes[1],"DNA")==0) {
			Maln.dnaflag=TRUE;
		}
	}
	else if(strcmp(Element,"aln-score")==0) {
		Maln.alnscore=atof(Content);
		Maln.validalnscore=TRUE;
	}
	else if(strcmp(Element,"seq-name")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].name,Content);
	}
	else if(strcmp(Element,"accession")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].access,Content);
	}
	else if(strcmp(Element,"nid")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].nid,Content);
	}
	else if(strcmp(Element,"definition")==0) {
		for(i=0,j=0;i<strlen(Content);i++)
			if(isalnum(Content[i]) || isspace(Content[i]) || Content[i]=='.' || Content[i]=='-') Maln.seqs[Fseq+Nseqs].title[j++]=Content[i];
	}
	else if(strcmp(Element,"organism")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].org,Content);
	}
	else if(strcmp(Element,"taxid")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].taxid,Content);
	}
	else if(strcmp(Element,"lifedomain")==0) {
		strcpy(Maln.seqs[Fseq+Nseqs].lifedomain,Content);
	}
	else if(strcmp(Element,"taxon")==0) {
		Ntaxons++;
		if(Ntaxons>=MAXTAXON) {
			warning("too many TAXON elements for sequence %s",Maln.seqs[Fseq+Nseqs].name);
		}
		else {
		l=strlen(Content);
		Maln.seqs[Fseq+Nseqs].taxon[Ntaxons]=(char *)ckalloc((l+1)*sizeof(char));
		strcpy(Maln.seqs[Fseq+Nseqs].taxon[Ntaxons],Content);
		}
	}
	else if(strcmp(Element,"surface-accessibility")==0) {
		if(Maln.seqs[Fseq+Nseqs].accessibility!=NULL) {
			ckfree(Maln.seqs[Fseq+Nseqs].accessibility);
			Maln.seqs[Fseq+Nseqs].accessibility=NULL;
		}
	}
	else if(strcmp(Element,"suracc-data")==0) {
		n=0;
		if(Content!=NULL) {
			for(i=0;i<strlen(Content);i++)
				if(isspace(Content[i])) n++;
			Maln.seqs[Fseq+Nseqs].accessibility=(float *)ckalloc((n+12)*sizeof(float));
			i=getfloatargs(Content,Maln.seqs[Fseq+Nseqs].accessibility,n+10);
			Maln.seqs[Fseq+Nseqs].naccessibility=i;
		}
	}
	else if(strcmp(Element,"residue-contact-list")==0) {
			NCresidue=0;
	}
	else if(strcmp(Element,"contact-residue1")==0) {
			Cresidue1=atoi(Content)-1;
			if(Rcontacts==NULL) {
				Rcontacts=(RCONTACT *)ckalloc((Cresidue1+10)*sizeof(RCONTACT));
				NC=Cresidue1;
			}
			else if(Cresidue1>NC) {
				Rcontacts=(RCONTACT *)ckrealloc(Rcontacts,(Cresidue1+10)*sizeof(RCONTACT));
				for(j=NC+1;j<Cresidue1+10;j++)
					Rcontacts[j].n=0;
				NC=Cresidue1;
			}
	}
	else if(strcmp(Element,"contact-residue2")==0) {
			NCresidue++;
			Rcontacts[Cresidue1].n=NCresidue;
			Rcontacts[Cresidue1].res[NCresidue-1]=atoi(Content)-1;
	}
	else if(strcmp(Element,"ec")==0) {
		strcpy(ec_text,Content);
                for(j=0;j<strlen(ec_text)-1;j++) {
                        if(!isdigit(ec_text[j])) ec_text[j]=' ';
                }
                n=sscanf(ec_text,"%d %d %d %d",&ec[0],&ec[1],&ec[2],&ec[3]);
                if(n!=4) {
                        for(l=0;l<4;l++) ec[l]=0;
                }
		for(l=0;l<4;l++) {
			Maln.seqs[Fseq+Nseqs].ec[l]=ec[l];
		}
	}
	else if(strcmp(Element,"subgroup")==0) {
		Maln.seqs[Fseq+Nseqs].simgroup=atoi(Content);
	}
	else if(strcmp(Element,"sense")==0) {
		Maln.seqs[Fseq+Nseqs].sense=atoi(Content);
	}
	else if(strcmp(Element,"hydrophobicity")==0) {
		Maln.seqs[Fseq+Nseqs].hydrophobicity=atof(Content);
	}
	else if(strcmp(Element,"fitem")==0) {
	}
	else if(strcmp(Element,"ftype")==0) {
		if(strcmp(Content,"DOMAIN")==0) 
			Ftype=SWDOMAIN;
		else if (strcmp(Content,"PFAM-A")==0) 
			Ftype=PFAMA;
		else if (strcmp(Content,"PFAM-B")==0) 
			Ftype=PFAMB;
		else if (strcmp(Content,"PROSITE")==0) 
			Ftype=PROSITE;
		else if (strcmp(Content,"ELM")==0) 
			Ftype=ELM;
		else if (strcmp(Content,"DISORDER")==0) 
			Ftype=DISORDER;
		else if (strcmp(Content,"SITE")==0) 
			Ftype=SITE;
		else if (strcmp(Content,"TRANSMEM")==0) 
			Ftype=TRANSMEM;
		else if (strcmp(Content,"REPEAT")==0) 
			Ftype=REPEAT;
		else if (strcmp(Content,"COIL")==0) 
			Ftype=COIL;
		else if (strcmp(Content,"SIGNAL")==0) 
			Ftype=SIGNAL;
		else if (strcmp(Content,"STRUCT")==0) 
			Ftype=STRUCT;
		else if (strcmp(Content,"VARSPLIC")==0) 
			Ftype=VARSPLIC;
		else if (strcmp(Content,"VARIANT")==0) 
			Ftype=VARIANT;
		else if (strcmp(Content,"ANCHOR")==0) 
			Ftype=ANCHOR;
		else if (strcmp(Content,"BLOCK")==0) 
			Ftype=COREBLOCK;
		else if (strcmp(Content,"SEQERR")==0) 
			Ftype=SEQERRBLOCK;
		else if (strcmp(Content,"REGION")==0) 
			Ftype=REGION;
		else if (strcmp(Content,"PHYLOBLOCK")==0) 
			Ftype=PHYLOBLOCK;
		else if (strcmp(Content,"MOD_RES")==0) 
			Ftype=MODRES;
		else if (strcmp(Content,"LOWCOMP")==0) 
			Ftype=LOWC;
		else if (strcmp(Content,"UNKNOWN")==0) 
			Ftype=OTHER;
		Nftable[Ftype]++;
		if(Nftable[Ftype]>=MAXFT) 
			warning("too many FTABLE elements for sequence %s",Maln.seqs[Fseq+Nseqs].name);
		else alloc_ft_entry(&Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]]);
		if(Nftable[Ftype]<MAXFT) strcpy(Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].type,Content);
	}
	else if(strcmp(Element,"fstart")==0) {
		if(Nftable[Ftype]<MAXFT) Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].start=atoi(Content)-1;
	}
	else if(strcmp(Element,"fstop")==0) {
		if(Nftable[Ftype]<MAXFT) Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].end=atoi(Content)-1;
	}
	else if(strcmp(Element,"fcolor")==0) {
		if(Nftable[Ftype]<MAXFT) Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].color=atoi(Content);
	}
	else if(strcmp(Element,"fscore")==0) {
		if(Nftable[Ftype]<MAXFT) Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].score=atof(Content);
	}
	else if(strcmp(Element,"fnote")==0) {
                if(Nftable[Ftype]<MAXFT)
                        if(Content!=NULL) {
                                for(i=0;i<strlen(Content);i++)
                                        if (!isspace(Content[i])) break;
                                strcpy(Maln.ft[Fseq+Nseqs].data[Ftype][Nftable[Ftype]].name,&Content[i]);
                        }
	}
	else if(strcmp(Element,"column-score")==0) {
		if(Ncolscores<MAXCSCORE) Ncolscores++;
	}
	else if(strcmp(Element,"colsco-name")==0) {
		if(Ncolscores<MAXCSCORE) if(Content!=NULL) {
			n=strlen(Content);
			Maln.col_score[Ncolscores].name=(char *)ckalloc((n+1)*sizeof(char));
			strcpy(Maln.col_score[Ncolscores].name,Content);
		}
	}
	else if(strcmp(Element,"colsco-owner")==0) {
		if(Ncolscores<MAXCSCORE) if(Content!=NULL) {
			n=strlen(Content);
			Maln.col_score[Ncolscores].owner=(char *)ckalloc((n+1)*sizeof(char));
			strcpy(Maln.col_score[Ncolscores].owner,Content);
		}
	}
	else if(strcmp(Element,"colsco-data")==0) {
		if(Ncolscores<MAXCSCORE) {
			n=0;
			if(Content!=NULL) {
				for(i=0;i<strlen(Content);i++)
					if(isspace(Content[i])) n++;
				Maln.col_score[Ncolscores].data=(sint *)ckalloc((n+12)*sizeof(sint));
				i=getintargs(Content,Maln.col_score[Ncolscores].data,n+10);
				Maln.col_score[Ncolscores].length=n;
				Maln.ncol_scores=Ncolscores+1;
			}
		}
	}
	else if(strcmp(Element,"goxref")==0) {
		Ngorefs++;
		if(Ngorefs>=MAXGOREF) 
			warning("too many GOXREF elements for sequence %s",Maln.seqs[Fseq+Nseqs].name);
	
	}
	else if(strcmp(Element,"goid")==0) {
		if(Ngorefs<MAXGOREF) {
		l=strlen(Content);
		Maln.go[Fseq+Nseqs].goref[Ngorefs].id=(void *)ckalloc((l+1)*sizeof(char));
		strcpy(Maln.go[Fseq+Nseqs].goref[Ngorefs].id,Content);
		}
	}
	else if(strcmp(Element,"goclass")==0) {
		if(Ngorefs<MAXGOREF) {
		Maln.go[Fseq+Nseqs].goref[Ngorefs].class=Content[0];
		}
	}
	else if(strcmp(Element,"godesc")==0) {
		if(Ngorefs<MAXGOREF) {
		l=strlen(Content);
		Maln.go[Fseq+Nseqs].goref[Ngorefs].desc=(void *)ckalloc((l+1)*sizeof(char));
		strcpy(Maln.go[Fseq+Nseqs].goref[Ngorefs].desc,Content);
		}
	}
	else if(strcmp(Element,"goevidence")==0) {
		if(Ngorefs<MAXGOREF) {
		strncpy(Maln.go[Fseq+Nseqs].goref[Ngorefs].evidence,Content,3);
		}
	}
	else if(strcmp(Element,"length")==0) {
		Maln.seqs[Fseq+Nseqs].len=atoi(Content);
	}
	else if(strcmp(Element,"group")==0) {
		Maln.seqs[Fseq+Nseqs].simgroup=atoi(Content);
	}
	else if(strcmp(Element,"keyword")==0) {
		Nkeywords++;
		if(Nkeywords>=MAXKEYWORDS) {
			warning("too many KEYWORD elements for sequence %s",Maln.seqs[Fseq+Nseqs].name);
		}
		else {
		l=strlen(Content);
		Maln.seqs[Fseq+Nseqs].keyword[Nkeywords]=(void *)ckalloc((l+1)*sizeof(char));
		strcpy(Maln.seqs[Fseq+Nseqs].keyword[Nkeywords],Content);
		}
	}
	else if(strcmp(Element,"seq-data")==0) {
		l=strlen(Content);
		alloc_seq(&Maln.seqs[Fseq+Nseqs],l);
		for(i=0,j=0;i<l;i++)
			if(!isspace(Content[i])) Maln.seqs[Fseq+Nseqs].data[j++]=Content[i];
		Maln.seqs[Fseq+Nseqs].len=j;
/* then this sequence must be finished, so we can update the data for this sequence */

		Maln.seqs[Fseq+Nseqs].nkeywords=Nkeywords+1;
		Maln.go[Fseq+Nseqs].ngorefs=Ngorefs+1;
		Maln.seqs[Fseq+Nseqs].ntaxons=Ntaxons+1;
		for(i=0;i<MAXFTTYPE;i++)
			Maln.ft[Fseq+Nseqs].nentries[i]=Nftable[i]+1;

/* copy residue contact data to multiple alignment taking into account gaps */
		if(Rcontacts!=NULL) {
			Maln.seqs[Fseq+Nseqs].rcontacts=(RCONTACT *)ckalloc((Maln.seqs[Fseq+Nseqs].len+1)*sizeof(RCONTACT));
			for(j=0;j<=Cresidue1;j++) {
				pos2col1(Maln.seqs[Fseq+Nseqs].data,j,&cpos1);
				Maln.seqs[Fseq+Nseqs].rcontacts[cpos1].n=Rcontacts[j].n;
				for(k=0;k<Rcontacts[j].n;k++) {
					pos2col1(Maln.seqs[Fseq+Nseqs].data,Rcontacts[j].res[k],&cpos2);
					Maln.seqs[Fseq+Nseqs].rcontacts[cpos1].res[k]=cpos2;
				}
			}
			ckfree(Rcontacts);
			Rcontacts=NULL;
		}
		Maln.nseqs++;
		Nseqs++;
	}

	for(i=0;i<NAttributes;i++) {
		ckfree(Attributes[i*2]);
		ckfree(Attributes[i*2+1]);
	}
	if(Element!=NULL) Element=ckfree(Element);
	if(Content!=NULL) Content=ckfree(Content);
}


/* This only gets called at the end of an element that contains data. Shouldn't need to touch it again! */
static void end(void *data, const char *el)
{
	if(Element!=NULL) {
		do_element();
	}
	Depth--;
}  

/* This reads the element names and puts them in Element. Shouldn't need to touch it again! */
static void start(void *data, const char *el, const char **attr)
{
	int i,n,len;

	if(Element!=NULL) {
		do_element();
	}

	len=strlen(el);
	Element=(char *)ckalloc((len+1)*sizeof(char));
	strcpy(Element,el);

	n=0;
	for (i = 0; attr[i]; i += 2) {
		n++;
	}
	NAttributes=n;

	if(NAttributes>0) {
		Attributes=(char **)ckalloc((NAttributes*2)*sizeof(char *));
		for (i = 0; attr[i]; i += 2) {
			len=strlen(attr[i]);
			Attributes[i]=(char *)ckalloc((len+1)*sizeof(char));
			strcpy(Attributes[i],attr[i]);
			len=strlen(attr[i+1]);
			Attributes[i+1]=(char *)ckalloc((len+1)*sizeof(char));
			strcpy(Attributes[i+1],attr[i+1]);
		}
	}
/*
	for (i = 0; attr[i]; i += 2) {
    printf(" %s='%s'", attr[i], attr[i + 1]);
	}
*/

	Depth++;
}  

/* This reads the contents of each element, and puts into Content. Shouldn't need to touch it again! */
static void charhndl(void *data, const char *text, int len)
{
	int i,l=0;

	if(len<=0) return;

	if(Content==NULL)
		Content=(char *)ckalloc((len+1)*sizeof(char));
	else {
		l=strlen(Content);
		Content=(char *)ckrealloc(Content,(l+len+2)*sizeof(char));
	}
	for (i = 0; i < len; i++)
		Content[l++] = text[i];
	if(Content[l-1]=='\n') Content[l-1]='\0';
	else Content[l]='\0';
}  

/* handlers used, if we're only counting the number of sequences */
static void c_start(void *data, const char *el, const char **attr)
{
	int i,len;

	if(strcmp(el,"sequence")==0) Nseqs++;

	Depth++;
}  

static void c_end(void *data, const char *el)
{
	Depth--;
}  

sint count_xml_seqs(FILE *fin)
{
  XML_Parser p;

  p = XML_ParserCreate(NULL);

  if (! p) {
    fprintf(stderr, "Couldn't allocate memory for parser\n");
    exit(-1);
  }

  XML_SetElementHandler(p, c_start, c_end);

  Nseqs=0;
  for (;;) {
    int done;
    int len;

    len = fread(Buff, 1, BUFFSIZE, fin);
    if (ferror(fin)) {
      fprintf(stderr, "Read error\n");
      exit(-1);
    }
    done = feof(fin);

    if (! XML_Parse(p, Buff, len, done)) {
      fprintf(stderr, "Parse error at line %d:\n%s\n",
	      XML_GetCurrentLineNumber(p),
	      XML_ErrorString(XML_GetErrorCode(p)));
      exit(-1);
    }

    if (done)
      break;
  }
  return Nseqs;
}

ALN read_xml(FILE *fin,int first_seq)
{
  XML_Parser p;

  p = XML_ParserCreate(NULL);

  if (! p) {
    fprintf(stderr, "Couldn't allocate memory for parser\n");
    exit(-1);
  }
  alloc_aln(Nseqs,&Maln);

  XML_SetElementHandler(p, start, end);
  XML_SetCharacterDataHandler(p, charhndl);

  Nseqs=0;
  Fseq=first_seq;
  for (;;) {
    int done;
    int len;

    len = fread(Buff, 1, BUFFSIZE, fin);
    if (ferror(fin)) {
      fprintf(stderr, "Read error\n");
      exit(-1);
    }
    done = feof(fin);

    if (! XML_Parse(p, Buff, len, done)) {
      fprintf(stderr, "Parse error at line %d:\n%s\n",
	      XML_GetCurrentLineNumber(p),
	      XML_ErrorString(XML_GetErrorCode(p)));
      exit(-1);
    }

    if (done)
      break;
  }
  return Maln;
}

int getfloatargs(char *inline1,float *args,int max)
{

        char    *inptr;
        char    *tstring;
/*
#ifndef MAC
        char    *strtok(char *s1, const char *s2);
#endif
*/
        int     i;

        inptr=inline1;
        for (i=0;i<=max;i++)
        {
                if ((tstring=strtok(inptr," \t\n"))==NULL)
                        break;
                args[i]=atof(tstring);
                inptr=NULL;
        }

        return(i);
}


