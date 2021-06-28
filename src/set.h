#ifndef HASH_SET_H
#define HASH_SET_H

#ifdef __cplusplus
extern "C" {
#endif

#define HASH_SET_SEED 26571997

//linked list node for hashing with chaining
typedef struct setnode_t {
  struct setnode_t* next;
  char* value;
  int sources;
} Setnode;

int set_insert(char* word, int length, int source, Setnode* set, int set_size);

int set_lookup(char* word, int length, Setnode* set, int set_size);

/**
 * Returns a char* array of all values stored in the set.
 */
char** get_values(Setnode *set, int set_size);

void set_deallocate(Setnode* set, int set_size);

void print_set(Setnode* set, int set_size);

int run_set_tests();


#ifdef __cplusplus
}
#endif

#endif
