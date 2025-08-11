
#include "allvars.h"
//#include "allvars_bt.h"



void determine_descendants(struct halo_catalogue *catA, struct halo_catalogue *catB, int entry, int snapnum);

int sort_candlist(const void *a, const void *b);
void prepare_index_list(struct halo_catalogue *cat);
void load_subhalo_catalogue(int num, struct halo_catalogue *cat);


void *mymalloc_bt(size_t n);
void myfree_bt(void *p);
void decide_upon_descendant(void);
void count_progenitors(struct halo_catalogue *catA, struct halo_catalogue *catB);

void save_decendant_list(void);

double second_bt(void);
double measure_time_bt(void);
double timediff_bt(double t0, double t1);

void generate_trees(void);

void load_subhalo_catalogue_ht(int num);
void count_halos(void);

void *mymalloc_ht(size_t n);
void myfree_ht(void *p);
void set_progenitor_pointers(void);

int peano_hilbert_key_ht(int x, int y, int z, int bits);
int whichfile(float *pos);


