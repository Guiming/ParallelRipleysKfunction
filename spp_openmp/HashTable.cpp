// Revised based on https://github.com/Karlina-Bytes/HashTable_Tutorial
// Revision made by Guiming Zhang (gzhang45@wisc.edu)

//*****************************************************************
//  HashTable.cpp
//  HashTable
//
//  Created by Kar Beringer on June 18, 2014.
//
//  This header file contains the Hash Table class definition.
//  Hash Table array elements consist of Linked List objects.
//*****************************************************************

#include "HashTable.h"
#include "Utility.h"
#include <cmath>

// Constructs the empty Hash Table object.
// Array length is set to 13 by default.
HashTable::HashTable(int tableLength, int grpFactor)
{
	if (tableLength <= 0) tableLength = 13;
	if (groupingFactor <= 0) groupingFactor = 10;
	array = new LinkedList[tableLength];
	length = tableLength;
	groupingFactor = grpFactor;

	// By Guiming, 2015-12-19
	keyLowerBound = 0.0;
	defaultItem = new Item(0.0, 1.0, NULL);
}

// Returns an array location for a given item key.
int HashTable::hash(double itemKey)
{
	return (int)(itemKey / groupingFactor) % length;
	//return (int)(itemKey) % length;
}

// Adds an item to the Hash Table if it is not there yet.
void HashTable::insertItem(Item * newItem)
{
	// By Guiming, 2015-12-19
	double w = newItem->weight;
	double d = newItem->distance;

	if (abs(w - 1.0) < EPSILON){ // no need to insert
		if (d > keyLowerBound) keyLowerBound = d;
		delete newItem;
		return;
	}

	int index = hash(d);
	if (array[index].getItem(index) == NULL){
		array[index].insertItem(newItem);
	}	
}

// Deletes an Item by key from the Hash Table.
// Returns true if the operation is successful.
bool HashTable::removeItem(double itemKey)
{
	// By Guiming, 2015-12-19
	if (itemKey <= keyLowerBound) return false;

	int index = hash(itemKey);
	return array[index].removeItem(itemKey);
}

// Deletes all Items from the Hash Table.
// Returns true if the operation is successful.
bool HashTable::clear()
{
	bool clearHashTab = true;
	bool clearLinkedList;
	for (int i = 0; i < length; i++){
		clearLinkedList = array[i].clear();
		clearHashTab = clearHashTab && clearLinkedList;
		if (!clearLinkedList)
			break;		
	}
	return clearHashTab;
}

// Returns an item from the Hash Table by key.
// If the item isn't found, a null pointer is returned.
Item * HashTable::getItemByKey(double itemKey)
{	
	// By Guiming, 2015-12-19
	if (itemKey <= keyLowerBound) {
		//Item item = Item(itemKey, 1.0, NULL);
		return defaultItem;
	}

	int index = hash(itemKey);
	return array[index].getItem(itemKey);
}

// Display the contents of the Hash Table to console window.
void HashTable::printTable()
{
	//cout << "Hash Table:\n";
	printf("Hash Table:\n");
	for (int i = 0; i < length; i++)
	{
		//cout << "Bucket " << i << ": ";
		printf("Bucket %d:", i);
		array[i].printList();
	}
}

// Prints a histogram illustrating the Item distribution.
void HashTable::printHistogram()
{
	//cout << "\nHash Table Contains ";
	printf("\nHash Table Contains ");
	//cout << getNumberOfItems() << " Items total.\n";
	printf("%d Items total.\n", getNumberOfItems());
	for (int i = 0; i < length; i++)
	{
		//cout << i + 1 << "\t";
		printf("%d\t", i + 1);
		for (int j = 0; j < array[i].getLength(); j++)
			//cout << " X";
			printf(" X");
		//cout << "\n";
		printf("\n");
	}
}

// Returns the number of locations in the Hash Table.
int HashTable::getLength()
{
	return length;
}

// Returns the number of Items in the Hash Table.
int HashTable::getNumberOfItems()
{
	int itemCount = 0;
	for (int i = 0; i < length; i++)
	{
		itemCount += array[i].getLength();
	}
	return itemCount;
}

// De-allocates all memory used for the Hash Table.
HashTable::~HashTable()
{
	delete[] array;
	delete defaultItem;
}

//*****************************************************************
// End of File
//*****************************************************************
