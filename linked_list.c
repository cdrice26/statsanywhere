// linked_list.c
#include "linked_list.h"

Node* create_linked_list(void) {
    Node *head = malloc(sizeof(Node));
    if (head == NULL) return NULL;
    head->data = NULL;
    head->next = NULL;
    return head;
}

int push(Node *node, void *data) {
    if (node == NULL) return 0;
    while (node->next != NULL) {
        node = node->next;
    }
    node->data = data;
    node->next = malloc(sizeof(Node));
    if (node->next == NULL) return 0;
    node->next->data = NULL;
    node->next->next = NULL;
    return 1;
}

// free_data is a callback for freeing whatever data points to,
// pass NULL if data doesn't need freeing (e.g. it's stack memory)
void free_linked_list(Node *head, void (*free_data)(void*)) {
    Node *temp = head;
    while (temp != NULL) {
        Node *next = temp->next;
        if (free_data && temp->data) free_data(temp->data);
        free(temp);
        temp = next;
    }
}

int length(Node *head) {
    int len = 0;
    Node *temp = head;
    while (temp->next != NULL) {
        len++;
        temp = temp->next;
    }
    return len;
}
