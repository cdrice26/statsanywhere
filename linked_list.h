#ifndef LINKED_LIST_H
#define LINKED_LIST_H
#include <stdlib.h>

// linked_list.h
typedef struct Node {
    void *data;
    struct Node *next;
} Node;

Node* create_linked_list(void);
int push(Node *head, void *data);
void free_linked_list(Node *head, void (*free_data)(void*));
int length(Node *head);

#endif
