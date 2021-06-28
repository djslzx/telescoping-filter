#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "set.h"
#include "murmur3.h"

//insert the string word into set
//return 1 if inserted in the set; 0 if already found (and not inserted)
int set_insert(char* word, int length, int source, Setnode* set, int set_size) {
  uint32_t hash[4];
  MurmurHash3_x64_128(word, length, HASH_SET_SEED, hash);
  int index = hash[0] % set_size;

  Setnode* currentNode = &set[index];

  if(currentNode->value == NULL) {
      currentNode->value = strdup(word);
      currentNode->sources = (1 << source);
  } else {
    while(currentNode->next != NULL){
      if(strcmp(currentNode->value, word) == 0) {
        currentNode->sources |= (1 << source);
        return 0;
      }
      currentNode = currentNode->next;
    }
    if(strcmp(currentNode->value, word) == 0) {
      currentNode->sources |= (1 << source);
      return 0;
    }
    Setnode* newNode = malloc(sizeof(*currentNode));
    currentNode->next = newNode;

    newNode->next = NULL;
    newNode->value = strdup(word);
    newNode->sources = (1 << source);
  }
  return 1;
}

//look up word in the set
//returns sources if it is found; returns 0 otherwise
int set_lookup(char* word, int length, Setnode* set, int set_size) {
  uint32_t hash[4];
  MurmurHash3_x64_128(word, length, HASH_SET_SEED, hash);
  int index = hash[0] % set_size;

  Setnode* currentNode = &set[index];
  if(currentNode->value == NULL) {
    return 0;
  } else {
    while(currentNode->next != NULL){
      if(strcmp(currentNode->value, word) == 0) {
        return currentNode->sources;
      }
      currentNode = currentNode->next;
    }
    if(strcmp(currentNode->value, word)){
      return 0;
    } else{
      return currentNode->sources;
    }
  }
}

void set_deallocate(Setnode* set, int set_size) {
  for(int i = 0; i < set_size; i++) {
    Setnode* current = &set[i];
    Setnode* next = current->next;
    while(next != NULL) {
      current = next;
      next = current->next;
      free(current->value);
      free(current);
    }
    free(set[i].value);
  }
  free(set);
}

char** get_values(Setnode *set, int set_size) {
  char** values = calloc(set_size, sizeof(char*));
  int j=0;
  for (int i=0; i<set_size; i++) {
    Setnode* node = &set[i];
    if (node->value != NULL) {
      do {
        values[j] = node->value;
        node = node->next;
        j++;
      } while (node != NULL);
    }
  }
  return values;
}

/* Printing */

void print_set(Setnode* set, int set_size) {
  printf("SET (size=%d):\n", set_size);
  for(int i = 0; i < set_size; i++) {
    if (set[i].value != NULL) {
      printf(" %d: [%s]", i, set[i].value);
      Setnode* next = set[i].next;
      while (next != NULL) {
        printf("-> [%s]", next->value);
        next = next->next;
      }
      printf("\n");
    }
  }
}

/* Tests */

/**
  * @return 1 if the set contains the word, 0 otherwise
 */
static int contains(char* set[], int set_size, char* word) {
  for (int i=0; i<set_size; i++)
    if (strcmp(set[i], word) == 0) return 1;
  return 0;
}

static int test_set_values() {
  printf("Starting get_values test...\n");
  char* inputs[] = {
          "premiere", "partie", "combray", "i", "longtemps", "me", "suis", "couche", "de", "bonne", "heure", "parfois",
          "a", "peine", "ma", "bougie", "eteinte", "mes", "yeux", "se", "fermaient", "si", "vite", "que", "n", "pas",
          "le", "temps", "dire", "m", "et", "une", "demi", "apres", "des", "epoques", "vecues", "par", "eux",
          "distantes", "entre", "lesquelles", "tant", "jours", "sont", "venus", "placer", "dans",
  };
  int num_words = sizeof(inputs) / sizeof(char*);
  int set_size = (int)(num_words * 1.5);
  Setnode* set = calloc(set_size, sizeof(set[0]));
  for (int i=0; i < num_words; i++) {
    set_insert(inputs[i], (int)strlen(inputs[i]), 0, set, set_size);
  }
  char** values = get_values(set, set_size);

  // Check that values = inputs by checking subsets both ways
  for (int i=0; i < num_words; i++) {
    if (!contains(values, num_words, inputs[i])) {
      printf("get_values test failed: input %s not found in set.\n", inputs[i]);
      return 0;
    }
    if (!contains(inputs, num_words, values[i])) {
      printf("get_values test failed: %s found in set but not in inputs.\n", values[i]);
      return 0;
    }
  }
  printf("get_values test passed.\n");
  return 1;
}

int run_set_tests() {
  return test_set_values();
}
