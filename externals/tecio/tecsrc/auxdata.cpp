#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE

/*
 *****************************************************************
 *****************************************************************
 *******                                                  ********
 ****** Copyright (C) 1988-2010 Tecplot, Inc.             ********
 *******       All Rights Reserved.                       ********
 *******                                                  ********
 *****************************************************************
 *****************************************************************
 */

#define AUXDATAMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "ALLOC.h"
#include "CHARTYPE.h"
#include "STRUTIL.h"
#include "ARRLIST.h"
#include "DATASET.h"
#include "STRLIST.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "SET.h"
#include "AUXDATA.h"

using namespace tecplot::strutil;

/**
 * Private auxiliary data item structure.
 */
typedef struct
{
    const char    *Name;
    ArbParam_t    Value;
    AuxDataType_e Type;
    Boolean_t     Retain;
} AuxDataItem_s;

/**
 * Private auxiliary data item container structure.
 */
typedef struct _AuxData_s
{
    /* invariant: ItemList is case-insensitive sorted by AuxDataItem->Name */
    ArrayList_pa ItemList; /* <AuxDataItem_s *>[dynamic] */
} AuxData_s;

//static Mutex_pa AuxDataMutex = NULL;


/**
 * A valid auxiliary data name character must begin with a '_' or alpha
 * character and may be followed by one or more '_', '.', alpha or digit
 * characters.
 */
Boolean_t AuxDataIsValidNameChar(char      Char,
                                 Boolean_t IsLeadChar)
{
    Boolean_t IsValidNameChar;

    REQUIRE("Char can be any value");
    REQUIRE(VALID_BOOLEAN(IsLeadChar));

    IsValidNameChar = (Char == '_' ||
                       tecplot::isalpha(Char));
    if (!IsLeadChar)
        IsValidNameChar = (IsValidNameChar ||
                           Char == '.'     ||
                           tecplot::isdigit(Char));

    ENSURE(VALID_BOOLEAN(IsValidNameChar));
    return IsValidNameChar;
}

/**
 * Indicates if the auxiliary data name is valid. A valid auxiliary data name
 * must begin with a '_' or alpha character and may be followed by one or
 * more '_', '.', alpha or digit characters.
 */
Boolean_t AuxDataIsValidName(const char *Name)
{
    Boolean_t  IsValidName;
    const char *NPtr;
    REQUIRE(VALID_REF(Name));

    for (NPtr = Name, IsValidName = AuxDataIsValidNameChar(*NPtr, TRUE);
         IsValidName && *NPtr != '\0';
         NPtr++)
    {
        IsValidName = AuxDataIsValidNameChar(*NPtr, FALSE);
    }

    ENSURE(VALID_BOOLEAN(IsValidName));
    return IsValidName;
}

/**
 * Deallocates an auxiliary data item and its contents and sets the target to
 * NULL.
 *
 * param AuxDataItem
 *     Reference to an auxiliary data item.
 */
static void AuxDataItemDealloc(AuxDataItem_s **AuxDataItem)
{
    REQUIRE(VALID_REF(AuxDataItem));
    REQUIRE(VALID_REF(*AuxDataItem) || *AuxDataItem == NULL);

    if (*AuxDataItem != NULL)
    {
        char *Name = (char *)(*AuxDataItem)->Name;
        if (Name != NULL)
            FREE_ARRAY(Name, "auxiliary name");

        if ((*AuxDataItem)->Type == AuxDataType_String)
        {
            char *Value = (char *)(*AuxDataItem)->Value;
            if (Value != NULL)
                FREE_ARRAY(Value, "auxiliary value");
        }
        else
            CHECK(FALSE);

        FREE_ITEM(*AuxDataItem, "auxiliary data item");
        *AuxDataItem = NULL;
    }

    ENSURE(*AuxDataItem == NULL);
}

/**
 * Allocates an auxiliary data item.
 *
 * NOTE: Copies are made of the name and value.
 *
 * param Name
 *     Auxiliary data item's name (case insenstive).
 * param Value
 *     Auxiliary data item's value.
 * param Type
 *     Auxiliary data item's value type.
 * param Retain
 *     Indicates if the auxiliary data item should persist. In other words
 *     copied, saved, etc.
 *
 * return
 *     A new auxiliary data item or NULL if sufficient memory was not
 *     available.
 */
static AuxDataItem_s *AuxDataItemAlloc(const char   *Name,
                                       ArbParam_t    Value,
                                       AuxDataType_e Type,
                                       Boolean_t     Retain)
{
    AuxDataItem_s *Result;

    REQUIRE(VALID_REF(Name) && AuxDataIsValidName(Name));
    REQUIRE(IMPLICATION(Type == AuxDataType_String,
                        (VALID_REF((char *)Value) ||
                         (char *)Value == NULL)));
    REQUIRE(VALID_ENUM(Type, AuxDataType_e));
    REQUIRE(VALID_BOOLEAN(Retain));

    Result = ALLOC_ITEM(AuxDataItem_s, "auxiliary data item");
    if (Result != NULL)
    {
        Boolean_t IsOk;
        Result->Type   = Type;
        Result->Retain = Retain;
        Result->Name   = DupString(dontTranslate(Name));
        IsOk = (Result->Name != NULL);
        Result->Value = 0; /* to satisfy some compilers' uninitialized warnings */
        if (IsOk && Type == AuxDataType_String)
        {
            char *strValue = (char *)Value;
            if (strValue != NULL)
            {
                char *strCopy = DupString(dontTranslate(strValue));
                Result->Value = (ArbParam_t)strCopy;
                IsOk = (strCopy != NULL);
            }
            else
                Result->Value = (ArbParam_t)NULL;
        }
        else
            CHECK(FALSE);

        if (!IsOk)
            AuxDataItemDealloc(&Result);
    }

    ENSURE(VALID_REF(Result) || Result == NULL);
    return Result;
}

/**
 * Destroys an auxiliary data item list item. This is an item destructor
 * callback for ArrayList's private data.
 *
 * param ItemRef
 *     Reference to the auxiliary data item to cleanup.
 * param ClientData
 *     Not used.
 *
 * return
 *     TRUE is a requirement
 */
static Boolean_t AuxDataItemListItemDestructor(void       *ItemRef,
                                               ArbParam_t  ClientData)
{
    AuxDataItem_s **AuxDataItemRef = (AuxDataItem_s **)ItemRef;

    REQUIRE(VALID_REF(AuxDataItemRef));
    REQUIRE(VALID_REF(*AuxDataItemRef) || *AuxDataItemRef == NULL);
    UNUSED(ClientData);

    if (*AuxDataItemRef != NULL)
        AuxDataItemDealloc(AuxDataItemRef);

    ENSURE(*AuxDataItemRef == NULL);
    return TRUE;
}

/**
 * Destroys an auxiliary data item. This is an item destructor
 * callback for ArrayList's private data.
 *
 * param ItemRef
 *     Reference to the auxiliary data to cleanup.
 * param ClientData
 *     Not used.
 *
 * return
 *     TRUE is a requirement
 */
Boolean_t AuxDataItemDestructor(void       *ItemRef,
                                ArbParam_t  ClientData)
{
    AuxData_pa *AuxDataRef = (AuxData_pa *)ItemRef;

    REQUIRE(VALID_REF(AuxDataRef));
    REQUIRE(VALID_REF(*AuxDataRef) || *AuxDataRef == NULL);
    UNUSED(ClientData);

    if (*AuxDataRef != NULL)
        AuxDataDealloc(AuxDataRef);

    ENSURE(*AuxDataRef == NULL);
    return TRUE;
}


/**
 * Duplicates an auxiliary data item if its Retain flag is TRUE or if directed
 * by the callback data. This is an item duplicator callback for ArrayList.
 *
 * param TargetItemRef
 *     Reference to the auxiliary data item to receive duplicate.
 * param SourceItemRef
 *     Reference to the auxiliary data item to duplicate.
 * param ClientData
 *     Boolean indicating if the Retain flag should be considered.
 *
 * return
 *     TRUE if the duplication was a success
 *     FALSE otherwise.
 */
static Boolean_t AuxDataItemDuplicator(void       *TargetItemRef,
                                       void       *SourceItemRef,
                                       ArbParam_t ClientData)
{
    Boolean_t IsOk = TRUE;
    AuxDataItem_s **TargetAuxDataItemRef = (AuxDataItem_s **)TargetItemRef;
    AuxDataItem_s **SourceAuxDataItemRef = (AuxDataItem_s **)SourceItemRef;
    Boolean_t       ConsiderRetain;

    REQUIRE(VALID_REF(TargetAuxDataItemRef));
    REQUIRE(VALID_REF(SourceAuxDataItemRef));
    REQUIRE(VALID_REF(*SourceAuxDataItemRef) || *SourceAuxDataItemRef == NULL);
    REQUIRE(VALID_BOOLEAN((Boolean_t)ClientData));

    ConsiderRetain = (Boolean_t)ClientData;

    /* duplicate the item */
    if (*SourceAuxDataItemRef != NULL &&
        (!ConsiderRetain || (*SourceAuxDataItemRef)->Retain))
    {
        *TargetAuxDataItemRef = AuxDataItemAlloc((*SourceAuxDataItemRef)->Name,
                                                 (*SourceAuxDataItemRef)->Value,
                                                 (*SourceAuxDataItemRef)->Type,
                                                 (*SourceAuxDataItemRef)->Retain);
        IsOk = (*TargetAuxDataItemRef != NULL);
    }
    else
        *TargetAuxDataItemRef = NULL;

    ENSURE(VALID_REF(*TargetAuxDataItemRef) || *TargetAuxDataItemRef == NULL);
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}

/**
 * Deallocates an auxiliary data handle and sets the handle to NULL.
 *
 * param AuxData
 *     Reference to an auxiliary data handle or reference to NULL.
 */
void LIBCALL AuxDataDealloc(AuxData_pa *AuxData)
{
    REQUIRE(VALID_REF(AuxData));
    REQUIRE(VALID_REF(*AuxData) || *AuxData == NULL);

    if (*AuxData != NULL)
    {
        ArrayListDealloc(&(*AuxData)->ItemList, AuxDataItemListItemDestructor, 0);
        FREE_ITEM(*AuxData, "auxiliary data container");
        *AuxData = NULL;
    }

    ENSURE(*AuxData == NULL);
}

/**
 * Allocates an auxiliary data handle.
 *
 * return
 *     Auxiliary data handle or NULL if sufficient memory was not available.
 */
AuxData_pa AuxDataAlloc(void)
{
    AuxData_pa Result = ALLOC_ITEM(AuxData_s, "auxiliary data container");
    if (Result != NULL)
    {
        Result->ItemList = ArrayListAlloc(0, ArrayListType_VoidPtr, NULL, 0);
        if (Result->ItemList == NULL)
            AuxDataDealloc(&Result);
    }

    ENSURE(VALID_REF(Result) || Result == NULL);
    return Result;
}

/**
 * Copies the auxiliary data and all its members who's Retain flag is TRUE
 * if the ConsiderRetain flag is TRUE otherwise it copies everything.
 */
AuxData_pa AuxDataCopy(AuxData_pa AuxData,
                       Boolean_t  ConsiderRetain)
{
    AuxData_pa Result;

    REQUIRE(VALID_REF(AuxData));
    REQUIRE(VALID_BOOLEAN(ConsiderRetain));

    Result = ALLOC_ITEM(AuxData_s, "auxiliary data container");
    if (Result != NULL)
    {
        Result->ItemList = ArrayListCopy(AuxData->ItemList,
                                         AuxDataItemDuplicator,
                                         ConsiderRetain);
        if (Result->ItemList != NULL)
        {
            if (ConsiderRetain)
            {
                /*
                 * Now pass through the array cleaning up the holes left by those
                 * auxiliary data item member who's Retain flag was FALSE and
                 * therefore left a VOID pointer because it was not copied.
                 */
                LgIndex_t ItemOffset = 0;
                LgIndex_t ItemCount = ArrayListGetCount(Result->ItemList);
                while (ItemOffset < ItemCount)
                {
                    /* if there is more than one in a row remove them all at once */
                    if (ArrayListGetVoidPtr(Result->ItemList, ItemOffset) == NULL)
                    {
                        LgIndex_t BaseOffsetToRemove = ItemOffset;
                        LgIndex_t NumItemsToRemove   = 1;
                        while (BaseOffsetToRemove + NumItemsToRemove < ItemCount &&
                               ArrayListGetVoidPtr(Result->ItemList,
                                                   BaseOffsetToRemove + NumItemsToRemove) == NULL)
                            NumItemsToRemove++;

                        /* delete the NULL items */
                        ArrayListDeleteItems(Result->ItemList,
                                             BaseOffsetToRemove,
                                             NumItemsToRemove,
                                             NULL, 0);

                        /*
                         * Update ItemCount but leave ItemOffset alone as it is now
                         * indexing the next item to examine.
                         */
                        ItemCount = ArrayListGetCount(Result->ItemList);
                    }
                    else
                        ItemOffset++;
                }
            }
        }
        else
            AuxDataDealloc(&Result);
    }

    ENSURE(VALID_REF(Result) || Result == NULL);
    return Result;
}

/**
 * Gets the current number of auxiliary data items maintained by the auxiliary.
 *
 * param AuxData
 *     Handle to auxiliary data.
 *
 * return
 *     Number of items maintained by the auxiliary data.
 */
LgIndex_t LIBCALL AuxDataGetNumItems(AuxData_pa AuxData)
{
    LgIndex_t NumItems;

    REQUIRE(VALID_REF(AuxData));

    NumItems = ArrayListGetCount(AuxData->ItemList);

    ENSURE(NumItems >= 0);
    return NumItems;
}

/**
 * Gets the item index of the name if found or if not found the index where an
 * auxiliary data item could be inserted.
 *
 * param AuxData
 *     Handle to auxiliary data.
 * param Name
 *     Name used for the search (case insensitive).
 * param ItemIndex
 *     Address to hold the index of the found item or the index where an
 *     auxiliary data item could be inserted.
 *
 * return
 *     TRUE if the named item was found,
 *     FALSE otherwise.
 */
Boolean_t AuxDataGetItemIndex(AuxData_pa  AuxData,
                              const char *Name,
                              LgIndex_t  *ItemIndex)
{
    Boolean_t FoundItem = FALSE;
    LgIndex_t Index;
    LgIndex_t NumItems;

    REQUIRE(VALID_REF(AuxData));
    INVARIANT("AuxData->ItemList is case-insensitive sorted by AuxDataItem->Name");
    REQUIRE(VALID_REF(Name) && AuxDataIsValidName(Name));
    REQUIRE(VALID_REF(ItemIndex));

    /*
     * Note that the current implementation just does a linear search
     * though the array looking for the index of the item or if not
     * found the index of the insertion point. This should be replaced
     * with a binary search.
     */
    NumItems = AuxDataGetNumItems(AuxData);

# if defined DO_LINEAR_SEARCH
    {
        for (Index = 0; Index < NumItems; Index++)
        {
            AuxDataItem_s *AuxDataItem =
                (AuxDataItem_s *)ArrayListGetVoidPtr(AuxData->ItemList, Index);
            int CompareResult = ustrcmp(AuxDataItem->Name, Name);
            if (CompareResult >= 0)
            {
                FoundItem = (CompareResult == 0);
                break;
            }
        }
    }
# else
    {
        int low, high;
        low = 0;
        high = NumItems - 1;
        Index = 0;
        while (low <= high)
        {
            AuxDataItem_s *AuxDataItem;
            int CompareResult;
            Index = (low + high) / 2;
            AuxDataItem = (AuxDataItem_s *)ArrayListGetVoidPtr(AuxData->ItemList, Index);
            CompareResult = ustrcmp(Name, AuxDataItem->Name);
            if (CompareResult < 0)
                high = Index - 1; /* If the new name is "less" than the one we're comparing to,
                                 don't change Index since Index is already in the right spot */
            else if (CompareResult > 0)
                low = ++Index; /* If the new name it "greater" than the one we're comparing
                              against, we want to make sure its Index is greater than
                              the current name's index as well, that's why we increment Index here. */
            else
            {
                FoundItem = TRUE;
                break;
            }
        }
    }
# endif

    *ItemIndex = Index;

    ENSURE(VALID_BOOLEAN(FoundItem));
    ENSURE(0 <= *ItemIndex &&
           ((FoundItem  && *ItemIndex <  NumItems) ||
            (!FoundItem && *ItemIndex <= NumItems)));
    return FoundItem;
}

/**
 * Gets the auxiliary data item at the specified index.
 *
 * NOTE: The name and value are a references, NOT copies.
 *
 * param AuxData
 *     Handle to auxiliary data.
 * param Index
 *     Index of the auxiliary data item of interest.
 * param Name
 *     Address to hold the auxiliary data item name.
 * param Value
 *     Address to hold the auxiliary data item value.
 * param Type
 *     Address to hold the auxiliary data item type.
 * param Retain
 *     Address to hold the auxiliary data item retain flag.
 */
void LIBCALL AuxDataGetItemByIndex(AuxData_pa    AuxData,
		                           LgIndex_t     Index,
		                           const char    **Name,
		                           ArbParam_t    *Value,
		                           AuxDataType_e *Type,
		                           Boolean_t     *Retain)
{
    AuxDataItem_s *AuxDataItem;

    REQUIRE(VALID_REF(AuxData));
    INVARIANT("AuxData->ItemList is case-insensitive sorted by AuxDataItem->Name");
    REQUIRE(0 <= Index && Index < ArrayListGetCount(AuxData->ItemList));
    REQUIRE(VALID_REF(Name));
    REQUIRE(VALID_REF(Value));
    REQUIRE(VALID_REF(Type));
    REQUIRE(VALID_REF(Retain));

    AuxDataItem = (AuxDataItem_s *)ArrayListGetVoidPtr(AuxData->ItemList, Index);
    *Name       = AuxDataItem->Name;
    *Value      = AuxDataItem->Value;
    *Type       = AuxDataItem->Type;
    *Retain     = AuxDataItem->Retain;

    ENSURE(VALID_REF(*Name) && AuxDataIsValidName(*Name));
    ENSURE(IMPLICATION(*Type == AuxDataType_String,
                       (VALID_REF((char *)(*Value)) ||
                        (char *)(*Value) == NULL)));
    ENSURE(VALID_ENUM(*Type, AuxDataType_e));
    ENSURE(VALID_BOOLEAN(*Retain));
}

/**
 * Gets the auxiliary data item by the specified name if it exists.
 *
 * NOTE: The name and value are a references, NOT copies.
 *
 * param AuxData
 *     Handle to auxiliary data.
 * param Name
 *     Name used for the search (case insensitive).
 * param Value
 *     Address to hold the auxiliary data item value.
 * param Type
 *     Address to hold the auxiliary data item type.
 * param Retain
 *     Address to hold the auxiliary data item retain flag.
 *
 * return
 *     TRUE if the an auxilary data item by the specified name was found,
 *     FALSE otherwise.
 */
Boolean_t AuxDataGetItemByName(AuxData_pa    AuxData,
                               const char    *Name,
                               ArbParam_t    *Value,
                               AuxDataType_e *Type,
                               Boolean_t     *Retain)
{
    Boolean_t FoundItem;
    LgIndex_t ItemIndex;

    REQUIRE(VALID_REF(AuxData));
    INVARIANT("AuxData->ItemList is case-insensitive sorted by AuxDataItem->Name");
    REQUIRE(VALID_REF(Name) && AuxDataIsValidName(Name));
    REQUIRE(VALID_REF(Value));
    REQUIRE(VALID_REF(Type));
    REQUIRE(VALID_REF(Retain));

    FoundItem = AuxDataGetItemIndex(AuxData, Name, &ItemIndex);
    if (FoundItem)
    {
        const char *SameName;
        AuxDataGetItemByIndex(AuxData, ItemIndex, &SameName,
                              Value, Type, Retain);
        CHECK(ustrcmp(Name, SameName) == 0);
    }

    ENSURE(VALID_BOOLEAN(FoundItem));
    ENSURE(IMPLICATION(FoundItem,
                       IMPLICATION(*Type == AuxDataType_String,
                                   (VALID_REF((char *)(*Value)) ||
                                    (char *)(*Value) == NULL))));
    ENSURE(IMPLICATION(FoundItem,
                       VALID_ENUM(*Type, AuxDataType_e)));
    ENSURE(IMPLICATION(FoundItem,
                       VALID_BOOLEAN(*Retain)));
    return FoundItem;
}


/**
 * Get a string value from AuxData and convert it to a boolean.
 */
Boolean_t AuxDataGetBooleanItemByName(AuxData_pa     AuxData, /* IN */
                                      const char    *Name,    /* IN */
                                      Boolean_t     *Value,   /* OUT */
                                      AuxDataType_e *Type,    /* OUT */
                                      Boolean_t     *Retain)  /* OUT */
{
    Boolean_t FoundItem;

    REQUIRE(VALID_REF(AuxData));
    INVARIANT("AuxData->ItemList is case-insensitive sorted by AuxDataItem->Name");
    REQUIRE(VALID_REF(Name) && AuxDataIsValidName(Name));
    REQUIRE(VALID_REF(Value));
    REQUIRE(VALID_REF(Type));
    REQUIRE(VALID_REF(Retain));

    ArbParam_t strValue;
    FoundItem = AuxDataGetItemByName(AuxData,
                                     Name,
                                     &strValue,
                                     Type,
                                     Retain);

    if (FoundItem &&
        (ustrcmp((char *)strValue, "YES")  == 0 ||
         ustrcmp((char *)strValue, "YEP")  == 0 ||
         ustrcmp((char *)strValue, "Y")    == 0 ||
         ustrcmp((char *)strValue, "TRUE") == 0 ||
         ustrcmp((char *)strValue, "T")    == 0 ||
         ustrcmp((char *)strValue, "ON")   == 0 ||
         ustrcmp((char *)strValue, "1")    == 0))
    {
        *Value = TRUE;
    }
    else
    {
        *Value = FALSE;
    }

    ENSURE(VALID_BOOLEAN(FoundItem));
    ENSURE(VALID_BOOLEAN(*Value));
    return FoundItem;
}


/**
 * Adds the auxiliary data item to the auxiliary data or replaces it if one
 * already exists by the same name.
 *
 * NOTE: The auxiliary data makes copies of the name and value.
 *
 * param AuxData
 *     Auxiliary data handle.
 * param Name
 *     Auxiliary data item's name (case insenstive).
 * param Value
 *     Auxiliary data item's value.
 * param Type
 *     Auxiliary data item's value type.
 * param Retain
 *     Indicates if the auxiliary data item should persist.
 *
 * return
 *     TRUE if the item was added to the auxiliary data.
 */
Boolean_t AuxDataSetItem(AuxData_pa    AuxData,
                         const char    *Name,
                         ArbParam_t    Value,
                         AuxDataType_e Type,
                         Boolean_t     Retain)
{
    Boolean_t     IsOk;
    AuxDataItem_s *AuxDataItem;

    REQUIRE(VALID_REF(AuxData));
    INVARIANT("AuxData->ItemList is case-insensitive sorted by AuxDataItem->Name");
    REQUIRE(VALID_REF(Name) && AuxDataIsValidName(Name));
    REQUIRE(IMPLICATION(Type == AuxDataType_String,
                        (VALID_REF((char *)Value) ||
                         (char *)Value == NULL)));
    REQUIRE(VALID_ENUM(Type, AuxDataType_e));
    REQUIRE(VALID_BOOLEAN(Retain));

    AuxDataItem = AuxDataItemAlloc(Name, Value, Type, Retain);
    IsOk = (AuxDataItem != NULL);
    if (IsOk)
    {
        LgIndex_t       ItemIndex;
        ArrayListItem_u ListItem;

        /* add or replace the item to the list */
        ListItem.VoidPtr = (void *)AuxDataItem;
        if (!AuxDataGetItemIndex(AuxData, Name, &ItemIndex))
            IsOk = ArrayListInsertItem(AuxData->ItemList, ItemIndex, ListItem);
        else
            IsOk = ArrayListSetItem(AuxData->ItemList, ItemIndex, ListItem,
                                    AuxDataItemListItemDestructor, 0);

        if (!IsOk)
            AuxDataItemDealloc(&AuxDataItem);
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    INVARIANT("AuxData->ItemList is case-insensitive sorted by AuxDataItem->Name");
    return IsOk;
}

/**
 * Deletes the auxiliary data item at the specified index.
 *
 * param AuxData
 *     Auxiliary data handle.
 * param Index
 *     Index of the auxiliary data item of interest.
 */
void AuxDataDeleteItemByIndex(AuxData_pa AuxData,
                              LgIndex_t  Index)
{
    REQUIRE(VALID_REF(AuxData));
    REQUIRE(0 <= Index && Index < ArrayListGetCount(AuxData->ItemList));

    ArrayListDeleteItem(AuxData->ItemList, Index, AuxDataItemListItemDestructor, 0);
}

/**
 * Deletes the auxiliary data item by the specified name if it exists.
 *
 * param AuxData
 *     Auxiliary data handle.
 * param Name
 *     Name used for the search (case insensitive).
 *
 * return
 *     TRUE if the an auxilary data item by the specified name was found,
 *     FALSE otherwise.
 */
Boolean_t AuxDataDeleteItemByName(AuxData_pa AuxData,
                                  const char *Name)
{
    Boolean_t FoundItem;
    LgIndex_t ItemIndex;

    REQUIRE(VALID_REF(AuxData));
    REQUIRE(VALID_REF(Name) && AuxDataIsValidName(Name));

    FoundItem = AuxDataGetItemIndex(AuxData, Name, &ItemIndex);
    if (FoundItem)
        AuxDataDeleteItemByIndex(AuxData, ItemIndex);

    ENSURE(VALID_BOOLEAN(FoundItem));
    return FoundItem;
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */
