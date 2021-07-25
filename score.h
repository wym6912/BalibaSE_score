/*#include "/us1/user/julie/dmalloc/malloc.h"*/
/*********************CLUSTALW.H*********************************************/
/****************************************************************************/

   /*
   Main header file for ClustalW.  Uncomment ONE of the following 4 lines
   depending on which compiler you wish to use.
   */

/*#define VMS 1                 VAX or ALPHA VMS */

/*#define MAC 1                 Think_C for Macintosh */

/*#define MSDOS 1               Turbo C for PC's */

#define UNIX 1                /*Ultrix/Decstation, Gnu C for 
                                Sun, IRIX/SGI, OSF1/ALPHA */

/***************************************************************************/
/***************************************************************************/


#include "general.h"

#define MAXNAMES		30	/* Max chars read for seq. names */
#define MAXTITLES		100      /* Title length */
#define MAXORGANISMS		100      /* Organism length */
#define MAXTAXID		20      /* Taxid length */
#define FILENAMELEN 	256             /* Max. file name length */
#define MAXFTNAMES              100     /* Max chars read for feature table names */
	
/* tree algorithms */
#define NJ	0
#define BIONJ	1
#define QUICKNJ	2

/* sequence file formats */
#define UNKNOWN   0
#define EMBLSWISS 1
#define PIR 	  2
#define PEARSON   3
#define GDE    	  4
#define CLUSTAL   5	/* DES */
#define MSF       6 /* DES */
#define RSF       7	/* JULIE */
#define USER      8	/* DES */
#define PHYLIP    9	/* DES */
#define NEXUS    10
#define GSCOPE	 11
#define RELACS	 12

/* secondary structure type in input file */
#define NONE      0
#define SECST     1
#define GMASK     2

/* used for profile 2 which may sometimes be just a list of unaligned sequences */
#define PROFILE 0
#define SEQUENCE 1

#define MAXOCOLORS 24
#define RED 0
#define BLUE 1
#define YELLOW 2
#define MAGENTA 4
#define VIOLET 13
#define ORANGE 5
#define TURQUOISE 18
#define LGREEN 19
#define GREEN 8
#define LRED 16 
#define VLGREEN 3
#define VLRED 9
#define DGRAY 24
#define LGRAY 25
#define GRAY 26

/* put bootstrap data in phylip output tree on tree nodes or branches (depends which tree drawing
   program you're using!) */
#define BS_NODE_LABELS 2
#define BS_BRANCH_LABELS 1

#define PAGE_LEN       22   /* Number of lines of help sent to screen */

#define PAGEWIDTH	80  /* maximum characters on output file page */
#define LINELENGTH     	60  /* Output file line length */
#define GCG_LINELENGTH 	50
#define GSCOPE_LINELENGTH 	100

#ifdef VMS						/* Defaults for VAX VMS */
#define COMMANDSEP '/'
#define DIRDELIM ']'		/* Last character before file name in full file specs */
#define INT_SCALE_FACTOR 1000 /* Scaling factor to convert float to integer for profile scores */

#elif MAC
#define COMMANDSEP '/'
#define DIRDELIM ':'
#define INT_SCALE_FACTOR 100  /* Scaling factor to convert float to integer for profile scores */

#elif MSDOS
#define COMMANDSEP '/'
#define DIRDELIM '\\'
#define INT_SCALE_FACTOR 100  /* Scaling factor to convert float to integer for profile scores */

#elif UNIX
#define COMMANDSEP '-'
#define DIRDELIM '/'
#define INT_SCALE_FACTOR 1000 /* Scaling factor to convert float to integer for profile scores */
#endif

#define NUMRES 26		/* max size of comparison matrix */
#define GAPCOL 30		/* position of gap open penalty in profile - must be >= NUMRES */
#define LENCOL 31		/* position of gap extension penalty in profile - must be >= NUMRES */
#define GAP1	'.'		/* code for gap introduced by alignment */
#define GAP2	'-'		/* code for gap read in with sequence */

#define INPUT 0			/* order of sequences in output file */
#define ALIGNED 1
#define TREEORDER 2

#define LEFT 1			/* direction in phylogenetic tree */
#define RIGHT 2

#define NODE 0			/* phylogenetic tree - LEAF is equivalent to a sequence */
#define LEAF 1

#define HELIX 1			/* used to look up score for helix on helix during alignment */
#define SHEET 2			/* used to look up score for sheet on sheet during alignment */

#define MOTIFBLK 1024		/* size of block for memory allocation for anchors */

#define NOW_GO_HOME 1

typedef struct {
        int color;
        int start;
        int end;
        int nseqs;
        int *seqs;
        char taxon[MAXORGANISMS];
        int euk;
        int bac;
        int arc;
        int vir;
} GROUP_NODE,*GROUP_NODEPTR;

typedef struct {
	char title[30];
	char string[30];
} MatMenuEntry;

typedef struct {
	int noptions;
	MatMenuEntry opt[10];
} MatMenu;

typedef struct {
        sint first;                     
        sint last;                     
	sint code;
} BLOCK,*BLOCKPTR;

typedef struct {
        sint **data;
	sint nseqs;
	sint len;
} PROF,*PROFPTR;

typedef struct {
	sint class;			/* ballast2 LMS class */
	sint chain;			/* ballast2 LMS chain */
	sint seq1;
	sint pos1;
	sint res1;
	sint seq2;
	sint pos2;
	sint res2;
	sint len;
	sint weight;
	char type; 			/* 0=motif in input file, 1=propagated motif */
} MOTIF,*MOTIF_PTR;

typedef struct {
	char filename[FILENAMELEN+1];	/* input sequence file name */
	char treename[FILENAMELEN+1];	/* input phylogenetic file name */
	sint nseqs;			/* number of sequences in the profile */
} PRF_DATA,*PRF_DATAPTR;

/*#define MAXKEYWORDS 50*/
#define MAXKEYWORDS 1
#define MAXTAXON 100

#define MAXRCONTACT 20

typedef struct {
	sint n;			/* number of residue contacts for a residue */
	sint res[MAXRCONTACT];	/* residue contacts  for a residue */
} RCONTACT,*RCONTACTPTR;

typedef struct {
	sint start;		/* index of first residue */
	sint len;		/* sequence length */
	sint reslen;		/* number of residues */
	sint weight;		/* sequence weight */
	sint simgroup;		/* groups by composition similarity */
	sint output_index;	/* position of the sequence in the output file */
	sint type;		/* 0 protein, 1 dna, 2 2d structure (prediction), 3 3d structure */
	float hydrophobicity;	/* average hydrophobicity score */
	Boolean fragment;	/* TRUE or FALSE, FALSE by default */
	sint ec[4];		/* EC number */
	char sense;		/* DNA sense */
	char *name;		/* sequence name */
	char *title;		/* comment line */
	char *org;		/* organism */
	char taxid[MAXTAXID];	/* organism taxid */
	char lifedomain[MAXORGANISMS];	/* life domain */
	sint ntaxons;		/* number of taxon */
	char *taxon[MAXTAXON];	/* lineage */
	char *access;		/* access number */
	char *nid;		/* id number */
	sint nkeywords;		/* number of keywords */
	char *keyword[MAXKEYWORDS];	/* keywords */
	float *accessibility;	/* accessibility score for each residue */
	sint naccessibility;	/* number of accessibility scores */
	RCONTACT *rcontacts;	/* list of residue contacts for each residue */
	char *data;		/* residues encoded as index into amino_acid_codes[] in readseq.c */
	char *mask;		/* a mask for anything you like */
	void *user_data;	/* anything you like! */
} SEQ,*SEQPTR;

#define KNOWN 0
#define PREDICTED 1
#define PROPAGATED 2

typedef struct {
	char *type;			/* type */
	sint start;			/* first residue */
	sint end;			/* last residue */
	char status;			/* status */
	sint start_col;			/* first residue */
	sint end_col;			/* last residue */
	sint color;			/* color */
	float score;			/* reliability score */
	char *name;			/* text string */
} FT_ENTRY,*FT_ENTRYPTR;

#define SWDOMAIN 0
#define PFAMA 1
#define PFAMB 2
#define REPEAT 3
#define PHYLOBLOCK 4
#define REGION 5
#define SEQERRBLOCK 6
#define VARSPLIC 7
#define LOWC 8
#define COREBLOCK 9
#define ANCHOR 10 
#define TRANSMEM 11
#define COIL 12
#define SIGNAL 13
#define STRUCT 14
#define PROSITE 15
#define SITE 16
#define MODRES 17 
#define VARIANT 18
#define ELM 19
#define DISORDER 20
#define OTHER 21
#define MAXFTTYPE 22
#define MAXCSCORE 100
#define MAXFT 1000

typedef struct {
	sint nentries[MAXFTTYPE];		/* number of entries in feature table for a sequence */
	FT_ENTRY data[MAXFTTYPE][MAXFT];		/* feature table data for a sequence */
} FT,*FTPTR;

typedef struct {
	char *name;			/* text string */
	char *owner;			/* text string */
	sint length;			/* length of score data */
	sint *data;			/* the column scores */
} CSCORE,*CSCOREPTR;


#define MAXGOREF 50

typedef struct {
	char *id;			/* id number */
	char class;			/* class (C,F,P) */
	char *desc;			/* description */
	char evidence[4];		/* evidence */
} GOREF,*GOREFPTR;

typedef struct {
	sint ngorefs;		/* number of GO database cross-references */
	GOREF goref[MAXGOREF];	/* GO database cross-references */
} GO,*GOPTR;

#define MAXREP 100

typedef struct {
	sint start;			/* first residue */
	sint end;			/* second residue */
	sint type;			/* length */
	sint offset;			/* weight */
	sint chain;			/* chain of lms's */
	sint index;			/* index to list of anchors */
} REP_ENTRY,*REP_ENTRYPTR;


typedef struct {
	sint nrepeats;			/* number of repeats for a sequence */
	REP_ENTRY data[MAXREP];		/* repeat data for a sequence */
} REP,*REPPTR;


typedef struct {
	sint len;			/* length of the group */
	sint *seqs;			/* sequences in the group */
} GROUP,*GROUPPTR;

typedef struct {
	sint ngroups;			/* number of sequence groups */
	GROUP *grp;			/* data for each group */
} SGROUP,*SGROUPPTR;

/* multiple alignment data (or sometimes profile data!) */
typedef struct {
	char filename[FILENAMELEN+1];	/* input sequence file name */
	char treename[FILENAMELEN+1];	/* input phylogenetic file name */
	char alphabet[NUMRES];		/* the alphabet */
	sint nseqs;			/* number of sequences */
	sint nanchors;			/* number of anchors */
	sint validalnscore;		/* 0 not calculated, 1 calculated */
	float alnscore;			/* score for the complete alignment */
	MOTIF *motifs;			/* anchors */
	Boolean dnaflag;		/* TRUE if sequences are DNA/RNA, FALSE if protein */
	SEQ  *seqs;			/* the sequences */
	PRF_DATA prf1;			/* if the sequences are actually organised as 2 profiles, put info. here */
	PRF_DATA prf2;
	FT  *ft;			/* the feature tables for each sequence */
	REP  *repeat;			/* repeat data for each sequence */
	GO   *go;			/* GO database cross-references */
	SGROUP groups;			/* groups by similarity */
	sint ncol_scores;		/* number of different column scores */
	CSCORE col_score[MAXCSCORE];	/* column scores */
	void *user_data;		/* anything you like! */
} ALN,*ALNPTR;

typedef struct {
	sint score[NUMRES][NUMRES];		/* residue comparison matrix */
} MATRIX,*MATRIXPTR;


/* command line and/or menu options, declared in param.h, initialised in interface.c */

/* secondary structure options */
typedef struct {
        sint            output_struct_penalties;
        Boolean         use_ss1;
        Boolean         use_ss2;
        sint        helix_penalty;
        sint        strand_penalty;
        sint        loop_penalty;
        sint        helix_end_minus;
        sint        helix_end_plus;
        sint        strand_end_minus;
        sint        strand_end_plus;
        sint        helix_end_penalty;
        sint        strand_end_penalty;
} SS_OPT,*SS_OPTPTR;

/* alignment output options */
typedef struct {
	Boolean         showaln;
	Boolean         output_clustal;
	Boolean         output_gcg;
	Boolean         output_phylip;
	Boolean         output_nbrf;
	Boolean         output_gde;
	Boolean         output_nexus;
	Boolean         output_gscope;
	Boolean         output_relacs;
	Boolean         output_tfa;
	Boolean         output_rsf;
	char		clustal_outname[FILENAMELEN+1];
	char		gcg_outname[FILENAMELEN+1];
	char		phylip_outname[FILENAMELEN+1];
	char		nbrf_outname[FILENAMELEN+1];
	char		gde_outname[FILENAMELEN+1];
	char		nexus_outname[FILENAMELEN+1];
	char		gscope_outname[FILENAMELEN+1];
	char		relacs_outname[FILENAMELEN+1];
	char		tfa_outname[FILENAMELEN+1];
	char		rsf_outname[FILENAMELEN+1];
	FILE		*clustal_outfile;
	FILE		*gcg_outfile;
	FILE		*phylip_outfile;
	FILE		*nbrf_outfile;
	FILE		*gde_outfile;
	FILE		*nexus_outfile;
	FILE		*gscope_outfile;
	FILE		*relacs_outfile;
	FILE		*tfa_outfile;
	FILE		*rsf_outfile;
	Boolean         lowercase; /* Flag for GDE output */
	Boolean         seq_numbers;
	sint		output_order;
	sint		output_names; /* 0 for same as input file, 1 for sequence id, 2 for sequence access */
} ALNOUT_OPT,*ALNOUT_OPTPTR;

/* quick pairwise options */
typedef struct {
	sint            dna_ktup;   /* parameters for DNA */
	sint            dna_wind_gap;
	sint            dna_signif;
	sint            dna_window;
	sint            prot_ktup;   /* parameters for proteins */
	sint            prot_wind_gap;
	sint            prot_signif;
	sint            prot_window;
	Boolean         percent;
} QUICKPW_OPT,*QUICKPW_OPTPTR;

/* long pairwise options */
typedef struct {
	float    dna_go_penalty;
	float    dna_ge_penalty;
	float    prot_go_penalty;
	float    prot_ge_penalty;
	float    transition_weight;
	sint     matnum;
	sint  	 treealgo;
	char     mtrxname[FILENAMELEN+1];
	sint     dnamatnum;
	char     dnamtrxname[FILENAMELEN+1];
	char     usermtrxname[FILENAMELEN+1];
	char     dnausermtrxname[FILENAMELEN+1];
	Boolean  add_motif;
	Boolean  quick_pairalign;
} PW_OPT,*PW_OPTPTR;

/* multiple alignment, protein gap options */
typedef struct {
	char            hyd_residues[20];
	sint            gap_dist;
	Boolean         no_hyd_penalties;
	Boolean         no_var_penalties;
	Boolean         no_pref_penalties;
	Boolean         use_endgaps;
	Boolean         nendgappenalties;
	Boolean         cendgappenalties;
} MULTGAP_OPT,*MULTGAP_OPTPTR;

/* multiple alignment options */
typedef struct {
	float           dna_gap_open;
	float           dna_gap_extend;
	float           transition_weight;
	sint            dnamatnum;
	char            dnamtrxname[FILENAMELEN+1];
	char            dnausermtrxname[FILENAMELEN+1];
	float           prot_gap_open;
	float           prot_gap_extend;
	sint            matnum;
	char            mtrxname[FILENAMELEN+1];
	char            usermtrxname[FILENAMELEN+1];
	sint            divergence_cutoff;
	Boolean         no_weights;
	Boolean         neg_matrix;
	Boolean         reset_alignments_new;
	Boolean         reset_alignments_all;
	Boolean		propogate_motifs;
	Boolean		use_motifs;
	Boolean		use_chains;
	Boolean		use_ss_motifs;
	char		motif_filename[FILENAMELEN+1];
	char		chain_filename[FILENAMELEN+1];
	MULTGAP_OPTPTR	gap_opt;
} MULT_OPT,*MULT_OPTPTR;

/* profile alignment options */
typedef struct {
	sint		profile_type;	/* PROFILE or SEQUENCE */
} PRF_OPT,*PRF_OPTPTR;

/* phylogenetic tree options */
typedef struct {
	sint  	 treealgo;
	Boolean         use_ambiguities;
	Boolean         tossgaps;
	Boolean         kimura;
	Boolean		use_matrix;
	char		matrixfile[FILENAMELEN];
	sint		cutoff;
	Boolean         output_tree_clustal;
	Boolean         output_tree_phylip;
	Boolean         output_tree_distances;
	Boolean         output_tree_nexus;
	char		clustal_outname[FILENAMELEN+1];
	char		phylip_outname[FILENAMELEN+1];
	char		dist_outname[FILENAMELEN+1];
} TREE_OPT,*TREE_OPTPTR;

/* bootstrap tree options */
typedef struct {
	sint            boot_ntrials;
	unsigned sint    boot_ran_seed;
	sint            bootstrap_format;
} BSTREE_OPT,*BSTREE_OPTPTR;

/* all the command line options */
typedef struct {
	sint explicit_type;
	SS_OPTPTR ss_opt;
	ALNOUT_OPTPTR alnout_opt;
	QUICKPW_OPTPTR quickpw_opt;
	PW_OPTPTR pw_opt;
	MULT_OPTPTR mult_opt;
	PRF_OPTPTR prf_opt;
	TREE_OPTPTR tree_opt;
	BSTREE_OPTPTR bstree_opt;
} OPT,*OPTPTR;

typedef struct {
	sint format;		/* 0=default blast format,1=bin */
	sint score[27][27];	/* matrix scores */
	sint gapscore[27];	/* gap scores in bin matrices */
	sint mat_avscore;	/* average mismatch score */
} COMP_MATRIX,*COMP_MATRIXPTR;

/* a series of matrices in a user file */
typedef struct {
	sint llimit;			/* lower limit of % identity between sequences */
	sint ulimit;			/* upper limit of % identity between sequences */
	char mat_filename[FILENAMELEN];	/* blast matrix filename */
} SeriesMat;

#define MAXMAT 10
typedef struct {
	Boolean user_series;		/* TRUE if the user file contains a series of matrices */
					/* FALSE if the user file is a single matrix in blast format */
	sint nmat;			/* number of matrices in the file */
	SeriesMat mat[MAXMAT];		/* the matrices */
} UserMatSeries;
	

typedef struct node {		/* phylogenetic tree structure used in readtree.c */
        struct node *left;
        struct node *right;
        struct node *parent;
        float dist;
        sint  leaf;
        int order;
        char name[64];
} INODE, *INODEPTR;


typedef struct {		/* phylogenetic tree structure used in readtree.c */
	Boolean rooted_tree;
	Boolean distance_tree;
	sint nleaves;
	sint nnodes;
	INODEPTR root;
	INODEPTR *leafptr;
} IN_TREE,*IN_TREEPTR;

/* A phylogenetic tree is described by a string "description", and the lengths
   "left_branch","right_branch".
   "description" contains 1 entry for each sequence,
	description[i]=1 	if sequence is a descendant of this node
		      =0	otherwise
*/

typedef struct{            /* phylogenetic tree structure used in trees.c */
        char *description;
        double  left_branch;
        double  right_branch;
} PTREE;

typedef struct {	/* count of number of residues in partial realignment block */
        sint n;		/* number of residues left on n-terminus */
        sint a;		/* number of residues in alignment block */
        sint c;		/* number of residues left on c-terminus */
} ALNCOUNT,*ALNCOUNTPTR;

/*
   Prototypes
*/

/* interface.c */
void init_options(OPTPTR opt);
sint parse_params(Boolean xmenus,char **args, OPTPTR opt,ALNPTR mult_aln,Boolean interactive);

/* readseq.c */
sint readseqs(char *filename,ALNPTR mult_aln);

/* util.c */

void alloc_aln(sint nseqs,ALNPTR mult_aln);
void realloc_aln(sint first_seq,sint nseqs,ALNPTR mult_aln);
void free_aln(ALNPTR mult_aln);
void alloc_seq(SEQ *seq,sint length);
void alloc_ft_entry(FT_ENTRY *ft_data);
void alloc_taxon_entry(char *taxon);
void alloc_go_entry(GOREF *goref);
void realloc_seq(SEQ *seq,sint length);
void free_seq(SEQ *seq);

void *ckalloc(size_t bytes);
void *ckrealloc(void *ptr, size_t bytes);
void *ckfree(void *ptr);
char prompt_for_yes_no(char *title,char *prompt);
void fatal(char *msg, ...);
void error(char *msg, ...);
void warning(char *msg, ...);
void info(char *msg, ...);
char *rtrim(char *str);
char *blank_to_(char *str);
char *upstr(char *str);
char *lowstr(char *str);
void getstr(char *instr, char *outstr);
double getreal(char *instr, double minx, double maxx, double def);
int getint(char *instr, int minx, int maxx, int def);
void do_system(void);
Boolean linetype(char *line, char *code);
Boolean keyword(char *line, char *code);
Boolean blankline(char *line);
void get_path(char *str, char *path);
int getargs(char *inline1, char *args[], int max);
int getintargs(char *inline1, sint *args, int max);
float countid(SEQ seq1,SEQ seq2);
float countid1(SEQ seq1,SEQ seq2);
void set_revision_level(char *revision_level);
char *get_revision_level(void);
void set_usemenu(Boolean flag);
Boolean get_usemenu(void);
void create_parameter_output(OPT opt,ALN mult_aln);
void check_fragments(ALNPTR mult_aln);
void remove_gap_pos(sint prf_no,ALNPTR mult_aln);
void sort_scores(float *scores,int f,int l);
void swap_scores(float *scores,int s1, int s2);
float normalise_score(float score,float n,float ntot,float ntotseq);
sint check_ft_type(char *ft_type,char *ft_name,sint *type);
void pos2col(char *seq,sint pstart,sint pend,sint *cstart,sint *cend);
void col2pos(char *seq,sint cstart,sint cend,sint *pstart,sint *pend);
void pos2col1(char *seq,sint pstart,sint *cstart);
void col2pos1(char *seq,sint cstart,sint *pstart);
sint overlap(sint f1,sint l1,sint f2, sint l2);
void add_ft_entry(ALNPTR mult_aln,sint seq,sint first,sint last,sint type,sint code,float score,char *ctype,char *name,sint is,sint ie);
void common_denominator(ALN mult_aln,GROUP_NODEPTR node);

/* readxml.c */
sint count_xml_seqs(FILE *fin);
ALN read_xml(FILE *fin,int first_seq);


