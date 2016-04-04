// Revised based on https://github.com/Karlina-Bytes/HashTable_Tutorial
// Revision made by Guiming Zhang (gzhang45@wisc.edu)

//*****************************************************************
//  LinkedList.cpp
//  HashTable
//
//  Created by Karlina Beringer on June 16, 2014.
//
//  This header file contains the Linked List class declaration.
//  Hash Table array elements consist of Linked List objects.
//*****************************************************************
#include <cmath>
#include "LinkedList.h"
#include "Utility.h"

//extern int IGNORABLEDIFF;	

// Constructs the empty linked list object.
// Creates the head node and sets length to zero.
LinkedList::LinkedList()
{
	head = new Item;
	head->next = NULL;
	length = 0;
}

// Inserts an item at the end of the list.
void LinkedList::insertItem(Item * newItem)
{
	if (!head->next)
	{
		head->next = newItem;
		length++;
		return;
	}

	/*
	Item * p = head;
	Item * q = head;
	while (q)
	{
		p = q;
		q = p->next;
	}
	p->next = newItem;
	newItem->next = NULL;
	length++;
	*/

	// should insert at the beginning
	newItem->next = head->next;
	head->next = newItem;
	length++;
}

// Removes an item from the list by item key.
// Returns true if the operation is successful.
bool LinkedList::removeItem(double itemKey)
{
	if (!head->next) return false;
	Item * p = head;
	Item * q = head;
	while (q)
	{
		//if (q->distance == itemKey)
		if (abs(q->distance - itemKey) < IGNORABLEDIFF)
		{
			p->next = q->next;
			delete q;
			length--;
			return true;
		}
		p = q;
		q = p->next;
	}
	return false;
}

// Removes an item from the list by item key.
// Returns true if the operation is successful.
bool LinkedList::clear()
{
	if (!head->next) return true;
	Item * p = head;
	while (!p)
	{
		head = head->next;
		delete p;
		p = head;
	}
	head = new Item;
	head->next = NULL;
	return true;
}

// Searches for an item by its key.
// Returns a reference to first match.
// Returns a NULL pointer if no match is found.
Item * LinkedList::getItem(double itemKey)
{	
	Item * p = head;
	Item * q = head;
	while (q)
	{
		p = q;
		if ((p != head) && abs(p->distance - itemKey) < IGNORABLEDIFF)
		//if (p != head)
			return p;
		q = p->next;
	}
	return NULL;
}

// Displays list contents to the console window.
void LinkedList::printList()
{
	if (length == 0)
	{
		//cout << "{ }\n";
		printf("{ }\n");
		return;
	}
	Item * p = head;
	Item * q = head;
	//cout << "{ ";
	printf("{ ");
	while (q)
	{
		p = q;
		if (p != head)
		{
			//cout << p->distance << ": " <<p->weight;
			printf("%2.1f: %2.1f",  p->distance, p->weight);
			if (p->next) 
				//cout << ", ";
				printf(", ");
			else 
				//cout << " ";
				printf(" ");
		}
		q = p->next;
	}
	//cout << "}\n";
	printf("}\n");
}

// Returns the length of the list.
int LinkedList::getLength()
{
	return length;
}

// De-allocates list memory when the program terminates.
LinkedList::~LinkedList()
{
	Item * p = head;
	Item * q = head;
	while (q)
	{
		p = q;
		q = p->next;
		if (q) delete p;
	}
}

//***************************************************************** 
// End of File 
//*****************************************************************
