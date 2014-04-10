#include "stdafx.h"
#include "MASTER.h"
/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/
#if !defined ARRLIST_h
#define ARRLIST_h

#if defined EXTERN
#  undef EXTERN
#endif
#if defined ARRLISTMODULE
#  define EXTERN
#else
#  define EXTERN extern
#endif

#if !defined TECPLOTKERNEL
typedef struct _ArrayList_s* ArrayList_pa;
#endif

typedef enum
{
    ArrayListType_UnsignedChar,
    ArrayListType_UnsignedShort,
    ArrayListType_UnsignedInt,
    ArrayListType_UnsignedLong,
    ArrayListType_Int64,
    ArrayListType_Char,
    ArrayListType_Short,
    ArrayListType_Int,
    ArrayListType_Long,
    ArrayListType_Float,
    ArrayListType_Double,
    ArrayListType_LgIndex,
    ArrayListType_EntIndex,
    ArrayListType_SmInteger,
    ArrayListType_Boolean,
    ArrayListType_ArbParam,
    ArrayListType_UnsignedCharPtr,
    ArrayListType_UnsignedShortPtr,
    ArrayListType_UnsignedIntPtr,
    ArrayListType_UnsignedLongPtr,
    ArrayListType_Int64Ptr,
    ArrayListType_CharPtr,
    ArrayListType_ShortPtr,
    ArrayListType_IntPtr,
    ArrayListType_LongPtr,
    ArrayListType_FloatPtr,
    ArrayListType_DoublePtr,
    ArrayListType_LgIndexPtr,
    ArrayListType_EntIndexPtr,
    ArrayListType_SmIntegerPtr,
    ArrayListType_BooleanPtr,
    ArrayListType_ArbParamPtr,
    ArrayListType_VoidPtr,
    ArrayListType_FunctionPtr,
    ArrayListType_Any,
    END_ArrayListType_e,
    ArrayListType_Invalid = BadEnumValue
} ArrayListType_e;

typedef union
{
    unsigned char   UnsignedChar;
    unsigned short  UnsignedShort;
    unsigned int    UnsignedInt;
    unsigned long   UnsignedLong;
    Int64_t         Int64;
    char            Char;
    short           Short;
    int             Int;
    long            Long;
    float           Float;
    double          Double;
    LgIndex_t       LgIndex;
    EntIndex_t      EntIndex;
    SmInteger_t     SmInteger;
    Boolean_t       BBoolean;  /* X-Windows uses Boolean */
    ArbParam_t      ArbParam;
    unsigned char*  UnsignedCharPtr;
    unsigned short* UnsignedShortPtr;
    unsigned int*   UnsignedIntPtr;
    unsigned long*  UnsignedLongPtr;
    Int64_t*        Int64Ptr;
    char*           CharPtr;
    short*          ShortPtr;
    int*            IntPtr;
    long*           LongPtr;
    float*          FloatPtr;
    double*         DoublePtr;
    LgIndex_t*      LgIndexPtr;
    EntIndex_t*     EntIndexPtr;
    SmInteger_t*    SmIntegerPtr;
    Boolean_t*      BooleanPtr;
    ArbParam_t*     ArbParamPtr;
    void*           VoidPtr;
    void (*FunctionPtr)(void);
} ArrayListItem_u;

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**
 * Visitor for traversing an array list. An iterator may not perform any
 * operation that will adjust the size of the list.  In other words it may not
 * insert or delete items from the list. However an iterator may perform a get
 * operation or a set operation that do not expand the list size.
 *
 * param ItemRef
 *     Reference to the array list item visited.
 * param ClientData
 *     Any client data required for the visitor.
 *
 * return
 *     TRUE to continue visiting items, otherwise
 *     FALSE to discontinue visiting
 */
typedef Boolean_t (*ArrayListItemVisitor_pf)(void*      ItemRef,
                                             ArbParam_t ClientData);
#if 0 /* use this stub as a starting place */
{
    REQUIRE(VALID_REF(TypeRef));
    REQUIRE(VALID_REF(*TypeRef) || *TypeRef == NULL);

    Boolean_t DoContinue = TRUE;
    <type>* TypeRef = static_cast<<type>*>(ItemRef);

    ENSURE(VALID_BOOLEAN(DoContinue));
    return DoContinue;
}
#endif


/**
 * Destructor for cleaning up one or more array list items. If a destructor is
 * not supplied then the array items are simply discarded.
 *
 * NOTE: The only change to ArrayListItemVisitor_pf is the policy which
 *       requires that the return value is TRUE.
 *
 * param ItemRef
 *     Reference to the array list item to cleanup.
 * param ClientData
 *     Any client data required for cleanup.
 *
 * return
 *     TRUE is a requirement
 */
typedef ArrayListItemVisitor_pf ArrayListItemDestructor_pf;


/**
 * Duplicator for copying one or more array list items. If a duplicator is not
 * supplied then the array items are simply copied. For pointer types this
 * means by reference. Note that if a duplicator is used the rules for
 * duplication and subsequent cleanup are defined by the client.
 *
 * param TargetItemRef
 *     Reference to the array list to receive the duplicate.
 * param SourceItemRef
 *     Reference to the array list item to duplicate.
 * param ClientData
 *     Any client data required for duplication.
 *
 * return
 *     TRUE if the duplication was a success
 *     FALSE otherwise. If the duplication failed it
 *     is the client's responsibility to cleanup any
 *     partial duplication
 */
typedef Boolean_t (*ArrayListItemDuplicator_pf)(void*      TargetItemRef,
                                                void*      SourceItemRef,
                                                ArbParam_t ClientData);
#if 0 /* use this stub as a starting place */
{
    REQUIRE(VALID_REF(TargetTypeRef));
    REQUIRE(VALID_REF(SourceTypeRef));
    REQUIRE(VALID_REF(*SourceTypeRef) || *SourceTypeRef == NULL);

    Boolean_t IsOk = TRUE;
    <type>* TargetTypeRef = static_cast<<type>*>(TargetItemRef);
    <type>* SourceTypeRef = static_cast<<type>*>(SourceItemRef);

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}
#endif


/**
 * Adjusts the capacity request as necessary to minimize memory reallocations
 * for large lists. Unless the request exceeds the maximum the adjusted
 * capacity will be at least as big as requested however it may be larger if it
 * is determined that the space requirement is growing faster. If the maximum
 * is exceeded zero should be returned.
 *
 * param ArrayList
 *     Array list requesting the change in capacity.
 * param CurrentCapacity
 *     Current capacity of the array list.
 * param RequestedCapacity
 *     Capacity request or zero for default size.
 * param ClientData
 *     Any client data needed for the adjustment.
 *
 * return
 *     Adjusted capacity that is at least as large as the request or zero if
 *     unable to satisfy the requested capacity.
 */
typedef LgIndex_t (*ArrayListCapacityRequestAdjuster_pf)(ArrayList_pa ArrayList,
                                                         LgIndex_t    CurrentCapacity,
                                                         LgIndex_t    RequestedCapacity,
                                                         ArbParam_t   ClientData);
#if 0 /* use this stub as a starting place */
{
    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE((RequestedCapacity == 0 && CurrentCapacity == 0) ||
            RequestedCapacity > ArrayList->Capacity);

    LgIndex_t Result;

    ENSURE(Result == 0 || Result >= RequestedCapacity);
    return Result;
}
#endif


/* private ArrayList structure: only exposed so STRUTIL can use it */
typedef struct _ArrayList_s
{
    char*            Array;           /* byte array for holding the items */
    ArrayListType_e  Type;            /* type of array items */
    SmInteger_t      ItemSize;        /* byte size of an individual item */
    LgIndex_t        Count;           /* number of items in the array */
    LgIndex_t        Capacity;        /* maximum holding capacity of the array */
    Boolean_t        IsVisitingItems; /* indicates if an iteration is in progress */
    ArrayListCapacityRequestAdjuster_pf CapacityRequestAdjuster;
    ArbParam_t                          CapacityRequestAdjusterClientData;
} ArrayList_s;


/**
 * Compares two array list elements. Note that either string may be
 * NULL as array lists allow for NULL elements.
 *
 * @param Item1
 *     Element to compare against Item2.
 * @param Item2
 *     Element to compare against Item1.
 * @param ClientData
 *     Contextual information that was passed to the 'ArrayListQSort' function.
 *
 * @return
 *     - A value less than zero if Item1 is less than Item2.
 *     - A value of zero if Item1 is equal to Item2.
 *     - A value greater than zero if Item1 is greater than Item2.
 */
typedef int (STDCALL *ArrayListItemComparator_pf)(ArrayListItem_u Item1,
                                                  ArrayListItem_u Item2,
                                                  ArbParam_t      ClientData);

EXTERN Boolean_t ArrayListIsValid(ArrayList_pa ArrayList);
EXTERN ArrayListType_e ArrayListGetType(ArrayList_pa ArrayList);
EXTERN Boolean_t ArrayListEnlargeCapacity(ArrayList_pa ArrayList,
                                          LgIndex_t    RequestedCapacity);
EXTERN ArrayList_pa ArrayListAlloc(LgIndex_t                           EstimatedCapacity,
                                   ArrayListType_e                     Type,
                                   ArrayListCapacityRequestAdjuster_pf CapacityRequestAdjuster = 0,
                                   ArbParam_t                          CapacityRequestAdjusterClientData = 0);
EXTERN void ArrayListDealloc(ArrayList_pa*              ArrayList,
                             ArrayListItemDestructor_pf ItemDestructor = 0,
                             ArbParam_t                 ClientData = 0);
EXTERN void ArrayListClear(ArrayList_pa ArrayList);
EXTERN void ArrayListDeleteAllItems(ArrayList_pa               ArrayList,
                                    ArrayListItemDestructor_pf ItemDestructor = 0,
                                    ArbParam_t                 ClientData = 0);
EXTERN void ArrayListDeleteItems(ArrayList_pa               ArrayList,
                                 LgIndex_t                  ItemOffset,
                                 LgIndex_t                  Count,
                                 ArrayListItemDestructor_pf ItemDestructor = 0,
                                 ArbParam_t                 ClientData = 0);
EXTERN void ArrayListDeleteItem(ArrayList_pa               ArrayList,
                                LgIndex_t                  ItemOffset,
                                ArrayListItemDestructor_pf ItemDestructor = 0,
                                ArbParam_t                 ClientData = 0);
EXTERN ArrayList_pa ArrayListRemoveItems(ArrayList_pa ArrayList,
                                         LgIndex_t    ItemOffset,
                                         LgIndex_t    Count);
EXTERN ArrayListItem_u ArrayListRemoveItem(ArrayList_pa ArrayList,
                                           LgIndex_t    ItemOffset);
EXTERN Boolean_t ArrayListInsertItem(ArrayList_pa    ArrayList,
                                     LgIndex_t       ItemOffset,
                                     ArrayListItem_u Item);
EXTERN Boolean_t ArrayListInsert(ArrayList_pa Target,
                                 LgIndex_t    ItemOffset,
                                 ArrayList_pa Source);
EXTERN Boolean_t ArrayListVisitItems(ArrayList_pa            ArrayList,
                                     LgIndex_t               ItemOffset,
                                     LgIndex_t               Count,
                                     ArrayListItemVisitor_pf ItemVisitor,
                                     ArbParam_t              ClientData);
EXTERN ArrayList_pa ArrayListGetItems(ArrayList_pa ArrayList,
                                      LgIndex_t    ItemOffset,
                                      LgIndex_t    Count);
EXTERN ArrayListItem_u ArrayListGetItem(ArrayList_pa ArrayList,
                                        LgIndex_t    ItemOffset);
EXTERN Boolean_t ArrayListSetItem(ArrayList_pa               ArrayList,
                                  LgIndex_t                  ItemOffset,
                                  ArrayListItem_u            Item,
                                  ArrayListItemDestructor_pf ItemDestructor = 0,
                                  ArbParam_t                 ClientData = 0);
EXTERN Boolean_t ArrayListAppendItem(ArrayList_pa    ArrayList,
                                     ArrayListItem_u Item);
EXTERN Boolean_t ArrayListAppend(ArrayList_pa Target,
                                 ArrayList_pa Source);
EXTERN ArrayList_pa ArrayListCopy(ArrayList_pa               ArrayList,
                                  ArrayListItemDuplicator_pf ItemDuplicator = 0,
                                  ArbParam_t                 ClientData = 0);
EXTERN void* ArrayListToArray(ArrayList_pa               ArrayList,
                              ArrayListItemDuplicator_pf ItemDuplicator,
                              ArbParam_t                 ClientData);
EXTERN ArrayList_pa ArrayListFromArray(void*                      Source,
                                       LgIndex_t                  Count,
                                       ArrayListType_e            Type,
                                       ArrayListItemDuplicator_pf ItemDuplicator = 0,
                                       ArbParam_t                 ClientData = 0);

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
EXTERN Boolean_t ArrayListBSearch(ArrayList_pa               ArrayList,
                                  ArrayListItem_u            Item,
                                  ArrayListItemComparator_pf Comparator,
                                  ArbParam_t                 ClientData,
                                  LgIndex_t*                 ItemIndex = 0);

#if defined USE_MACROS_FOR_FUNCTIONS
/**
 * Gets the array list's internal buffer representation.
 * Use ArrayListGetXxx accessors whenever possible as their
 * implementation in the release build is as fast as using
 * the array directly anyway.
 *
 * WARNING:
 *     Some array list functions modify the internal buffer.
 *     This will invalidate this reference however it is
 *     the client's responsibility not to make further use
 *     of it. In addition, this reference should never be
 *     deallocated directly as the array list assumes the
 *     responsible for the cleanup.
 *
 * param ArrayList
 *     Array list for which a reference to the internal
 *     buffer is desired.
 *
 * return
 *     Reference to the array list's internal buffer.
 */
#  define ArrayListGetInternalRef     ArrayListGetInternalRef_MACRO
/**
 * Gets the item's internal reference at the specified offset in the list.
 *
 * WARNING:
 *     Some array list functions modify the internal buffer.
 *     This will invalidate this reference however it is
 *     the client's responsibility not to make further use
 *     of it. In addition, this reference should never be
 *     deallocated directly as the array list assumes the
 *     responsible for the cleanup.
 *
 * param ArrayList
 *     Array list containing the desired item.
 * param ItemOffset
 *     Offset to the item in the list.
 *
 * return
 *     The internal reference to the requested item.
 */
#  define ArrayListGetItemInternalRef ArrayListGetItemInternalRef_MACRO
#  define ArrayListGetCount           ArrayListGetCount_MACRO

#  define ArrayListGetUnsignedChar(ArrayList, ItemOffset)     ArrayListGetTypedItem(ArrayList, ItemOffset, unsigned char)
#  define ArrayListGetUnsignedShort(ArrayList, ItemOffset)    ArrayListGetTypedItem(ArrayList, ItemOffset, unsigned short)
#  define ArrayListGetUnsignedInt(ArrayList, ItemOffset)      ArrayListGetTypedItem(ArrayList, ItemOffset, unsigned int)
#  define ArrayListGetUnsignedLong(ArrayList, ItemOffset)     ArrayListGetTypedItem(ArrayList, ItemOffset, unsigned long)
#  define ArrayListGetInt64(ArrayList, ItemOffset)            ArrayListGetTypedItem(ArrayList, ItemOffset, Int64_t)
#  define ArrayListGetChar(ArrayList, ItemOffset)             ArrayListGetTypedItem(ArrayList, ItemOffset, char)
#  define ArrayListGetShort(ArrayList, ItemOffset)            ArrayListGetTypedItem(ArrayList, ItemOffset, short)
#  define ArrayListGetInt(ArrayList, ItemOffset)              ArrayListGetTypedItem(ArrayList, ItemOffset, int)
#  define ArrayListGetLong(ArrayList, ItemOffset)             ArrayListGetTypedItem(ArrayList, ItemOffset, long)
#  define ArrayListGetFloat(ArrayList, ItemOffset)            ArrayListGetTypedItem(ArrayList, ItemOffset, float)
#  define ArrayListGetDouble(ArrayList, ItemOffset)           ArrayListGetTypedItem(ArrayList, ItemOffset, double)
#  define ArrayListGetLgIndex(ArrayList, ItemOffset)          ArrayListGetTypedItem(ArrayList, ItemOffset, LgIndex_t)
#  define ArrayListGetEntIndex(ArrayList, ItemOffset)         ArrayListGetTypedItem(ArrayList, ItemOffset, EntIndex_t)
#  define ArrayListGetSmInteger(ArrayList, ItemOffset)        ArrayListGetTypedItem(ArrayList, ItemOffset, SmInteger_t)
#  define ArrayListGetBoolean(ArrayList, ItemOffset)          ArrayListGetTypedItem(ArrayList, ItemOffset, Boolean_t)
#  define ArrayListGetArbParam(ArrayList, ItemOffset)         ArrayListGetTypedItem(ArrayList, ItemOffset, ArbParam_t)
#  define ArrayListGetUnsignedCharPtr(ArrayList, ItemOffset)  ArrayListGetTypedItem(ArrayList, ItemOffset, unsigned char*)
#  define ArrayListGetUnsignedShortPtr(ArrayList, ItemOffset) ArrayListGetTypedItem(ArrayList, ItemOffset, unsigned short*)
#  define ArrayListGetUnsignedIntPtr(ArrayList, ItemOffset)   ArrayListGetTypedItem(ArrayList, ItemOffset, unsigned int*)
#  define ArrayListGetUnsignedLongPtr(ArrayList, ItemOffset)  ArrayListGetTypedItem(ArrayList, ItemOffset, unsigned long*)
#  define ArrayListGetInt64Ptr(ArrayList, ItemOffset)         ArrayListGetTypedItem(ArrayList, ItemOffset, Int64_t*)
#  define ArrayListGetCharPtr(ArrayList, ItemOffset)          ArrayListGetTypedItem(ArrayList, ItemOffset, char*)
#  define ArrayListGetShortPtr(ArrayList, ItemOffset)         ArrayListGetTypedItem(ArrayList, ItemOffset, short*)
#  define ArrayListGetIntPtr(ArrayList, ItemOffset)           ArrayListGetTypedItem(ArrayList, ItemOffset, int*)
#  define ArrayListGetLongPtr(ArrayList, ItemOffset)          ArrayListGetTypedItem(ArrayList, ItemOffset, long*)
#  define ArrayListGetFloatPtr(ArrayList, ItemOffset)         ArrayListGetTypedItem(ArrayList, ItemOffset, float*)
#  define ArrayListGetDoublePtr(ArrayList, ItemOffset)        ArrayListGetTypedItem(ArrayList, ItemOffset, double*)
#  define ArrayListGetLgIndexPtr(ArrayList, ItemOffset)       ArrayListGetTypedItem(ArrayList, ItemOffset, LgIndex_t*)
#  define ArrayListGetEntIndexPtr(ArrayList, ItemOffset)      ArrayListGetTypedItem(ArrayList, ItemOffset, EntIndex_t*)
#  define ArrayListGetSmIntegerPtr(ArrayList, ItemOffset)     ArrayListGetTypedItem(ArrayList, ItemOffset, SmInteger_t*)
#  define ArrayListGetBooleanPtr(ArrayList, ItemOffset)       ArrayListGetTypedItem(ArrayList, ItemOffset, Boolean_t*)
#  define ArrayListGetArbParamPtr(ArrayList, ItemOffset)      ArrayListGetTypedItem(ArrayList, ItemOffset, ArbParam_t*)
#  define ArrayListGetVoidPtr(ArrayList, ItemOffset)          ArrayListGetTypedItem(ArrayList, ItemOffset, void*)
#  define ArrayListGetFunctionPtr(ArrayList, ItemOffset)      ArrayListGetTypedItem(ArrayList, ItemOffset, (**)(void))
#else
#  define ArrayListGetInternalRef     ArrayListGetInternalRef_FUNC
#  define ArrayListGetItemInternalRef ArrayListGetItemInternalRef_FUNC
#  define ArrayListGetCount           ArrayListGetCount_FUNC

#  define ArrayListGetUnsignedChar(ArrayList, ItemOffset)       (*(static_cast<unsigned char*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetUnsignedShort(ArrayList, ItemOffset)     (*(static_cast<unsigned short*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetUnsignedInt(ArrayList, ItemOffset)         (*(static_cast<unsigned int*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetUnsignedLong(ArrayList, ItemOffset)       (*(static_cast<unsigned long*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetInt64(ArrayList, ItemOffset)                    (*(static_cast<Int64_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetChar(ArrayList, ItemOffset)                        (*(static_cast<char*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetShort(ArrayList, ItemOffset)                      (*(static_cast<short*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetInt(ArrayList, ItemOffset)                          (*(static_cast<int*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetLong(ArrayList, ItemOffset)                        (*(static_cast<long*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetFloat(ArrayList, ItemOffset)                      (*(static_cast<float*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetDouble(ArrayList, ItemOffset)                    (*(static_cast<double*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetLgIndex(ArrayList, ItemOffset)                (*(static_cast<LgIndex_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetEntIndex(ArrayList, ItemOffset)              (*(static_cast<EntIndex_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetSmInteger(ArrayList, ItemOffset)            (*(static_cast<SmInteger_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetBoolean(ArrayList, ItemOffset)                (*(static_cast<Boolean_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetArbParam(ArrayList, ItemOffset)              (*(static_cast<ArbParam_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetUnsignedCharPtr(ArrayList, ItemOffset)   (*(static_cast<unsigned char**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetUnsignedShortPtr(ArrayList, ItemOffset) (*(static_cast<unsigned short**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetUnsignedIntPtr(ArrayList, ItemOffset)     (*(static_cast<unsigned int**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetUnsignedLongPtr(ArrayList, ItemOffset)   (*(static_cast<unsigned long**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetInt64Ptr(ArrayList, ItemOffset)                (*(static_cast<Int64_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetCharPtr(ArrayList, ItemOffset)                    (*(static_cast<char**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetShortPtr(ArrayList, ItemOffset)                  (*(static_cast<short**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetIntPtr(ArrayList, ItemOffset)                      (*(static_cast<int**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetLongPtr(ArrayList, ItemOffset)                    (*(static_cast<long**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetFloatPtr(ArrayList, ItemOffset)                  (*(static_cast<float**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetDoublePtr(ArrayList, ItemOffset)                (*(static_cast<double**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetLgIndexPtr(ArrayList, ItemOffset)            (*(static_cast<LgIndex_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetEntIndexPtr(ArrayList, ItemOffset)          (*(static_cast<EntIndex_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetSmIntegerPtr(ArrayList, ItemOffset)        (*(static_cast<SmInteger_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetBooleanPtr(ArrayList, ItemOffset)            (*(static_cast<Boolean_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetArbParamPtr(ArrayList, ItemOffset)          (*(static_cast<ArbParam_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetVoidPtr(ArrayList, ItemOffset)                    (*(static_cast<void**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#  define ArrayListGetFunctionPtr(ArrayList, ItemOffset)             (*(static_cast<**(void)*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(ArrayList, ItemOffset)))))
#endif

#if !defined USE_MACROS_FOR_FUNCTIONS
EXTERN void const* ArrayListGetInternalRef_FUNC(ArrayList_pa ArrayList);
EXTERN void const* ArrayListGetItemInternalRef_FUNC(ArrayList_pa ArrayList,
                                                    LgIndex_t    ItemOffset);
EXTERN LgIndex_t ArrayListGetCount_FUNC(ArrayList_pa ArrayList);
#endif

#define ArrayListGetInternalRef_MACRO(ArrayList)                 static_cast<void const*>((ArrayList)->Array)
#define ArrayListGetItemInternalRef_MACRO(ArrayList, ItemOffset) static_cast<void const*>(&((ArrayList)->Array[(ItemOffset)*(ArrayList)->ItemSize]))
#define ArrayListGetCount_MACRO(ArrayList)                       ((ArrayList)->Count)
#define ArrayListGetTypedArrayRef(ArrayList, NativeType)         reinterpret_cast<NativeType*>((ArrayList)->Array)
#define ArrayListGetTypedItem(ArrayList, ItemOffset, NativeType) (ArrayListGetTypedArrayRef(ArrayList,NativeType)[ItemOffset])

/**
 * Simple macro to determine if the item offset is within the array list
 * capacity. In the debug or checked builds we also perform a lower bound
 * assertion check.
 */
#if defined NO_ASSERTS
# define ArrayListOffsetWithinCapacity(ArrayList, ItemOffset) ((ItemOffset) < (ArrayList)->Capacity)
#else
/**
 * Using 'assert' rather than 'REQUIRE' because under Windows, REQUIRE (and ASSERT) trickles down to being a
 * do-while loop, which doesn't jive well with the comma operator.
 */
# define ArrayListOffsetWithinCapacity(ArrayList, ItemOffset) ((assert((ItemOffset) >= 0),TRUE) && ((ItemOffset) < (ArrayList)->Capacity))
#endif

/**
 * Places the item at the specified offset. If the offset is beyond the
 * end of the list it is sized accordingly and the intervening items
 * between the last item of the original state and the last item of the
 * new state are guaranteed to be 0.
 *
 * This is the workhorse of the set and append convenience macros that follow.
 * Please note that unlike ArrayListSetItem no destructor facility is provided
 * therefore if an item previously occupied 'ItemOffset' it will be replaced.
 *
 * param ArrayList
 *     Array list target in which to set the item.
 * param ItemOffset
 *     Offset of the item.
 * param Item
 *     Item to set at the specified offset. Its native type must
 *     match 'NativeType'
 * param NativeType
 *     Native type of 'Item'.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
#define ArrayListSetTypedItem(ArrayList, ItemOffset, Item, NativeType) \
    ((ArrayListOffsetWithinCapacity((ArrayList), (ItemOffset)) || \
      ArrayListEnlargeCapacity((ArrayList), (ItemOffset)+1)) \
          ? ((void)((ArrayListGetTypedArrayRef((ArrayList),NativeType)[(ItemOffset)]) = (Item)), \
             (((ItemOffset)+1 > (ArrayList)->Count) \
                  ? (((ArrayList)->Count = (ItemOffset)+1), TRUE) \
                  : (TRUE))) \
          : (FALSE))

/**
 * Appends the item to the end of the list. This is similar to ArrayListSetTypedItem however it
 * needs not perform a bounds check which is always one beyond the end. This macro was added
 * primarily to remove compiler warnings caused by a comparison of the array list ArrayList->Count+1
 * always returning true when tested if larger than ArrayList->Count which occurs if you use the
 * ArrayListSetTypedItem as ArrayListSetTypedItem(ArrayList, (ArrayList)->Count, Item, NativeType).
 *
 * This is the workhorse of the append convenience macros that follow.
 *
 * param ArrayList
 *     Array list target in which to set the item.
 * param Item
 *     Item to append to the end of the array. Its native type must
 *     match 'NativeType'
 * param NativeType
 *     Native type of 'Item'.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
#define ArrayListAppendTypedItem(ArrayList, Item, NativeType) \
    ((ArrayListOffsetWithinCapacity((ArrayList), (ArrayList)->Count) || \
      ArrayListEnlargeCapacity((ArrayList), (ArrayList)->Count+1)) \
          ? ((void)((ArrayListGetTypedArrayRef((ArrayList),NativeType)[(ArrayList)->Count]) = (Item)), \
             (((ArrayList)->Count = (ArrayList)->Count+1), TRUE)) \
          : (FALSE))

/**
 * This section provides macros for high speed setting and appending to an
 * array list of a known type. The only additional overhead incurred versus just
 * using a simple array is the cost of testing the array list capacity.
 */
#define ArrayListSetUnsignedChar(ArrayList, ItemOffset, Item)     ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),unsigned char)
#define ArrayListSetUnsignedShort(ArrayList, ItemOffset, Item)    ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),unsigned short)
#define ArrayListSetUnsignedInt(ArrayList, ItemOffset, Item)      ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),unsigned int)
#define ArrayListSetUnsignedLong(ArrayList, ItemOffset, Item)     ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),unsigned long)
#define ArrayListSetInt64(ArrayList, ItemOffset, Item)            ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),Int64_t)
#define ArrayListSetChar(ArrayList, ItemOffset, Item)             ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),char)
#define ArrayListSetShort(ArrayList, ItemOffset, Item)            ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),short)
#define ArrayListSetInt(ArrayList, ItemOffset, Item)              ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),int)
#define ArrayListSetLong(ArrayList, ItemOffset, Item)             ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),long)
#define ArrayListSetFloat(ArrayList, ItemOffset, Item)            ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),float)
#define ArrayListSetDouble(ArrayList, ItemOffset, Item)           ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),double)
#define ArrayListSetLgIndex(ArrayList, ItemOffset, Item)          ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),LgIndex_t)
#define ArrayListSetEntIndex(ArrayList, ItemOffset, Item)         ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),EntIndex_t)
#define ArrayListSetSmInteger(ArrayList, ItemOffset, Item)        ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),SmInteger_t)
#define ArrayListSetBoolean(ArrayList, ItemOffset, Item)          ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),Boolean_t)
#define ArrayListSetArbParam(ArrayList, ItemOffset, Item)         ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),ArbParam_t)
#define ArrayListSetUnsignedCharPtr(ArrayList, ItemOffset, Item)  ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),unsigned char*)
#define ArrayListSetUnsignedShortPtr(ArrayList, ItemOffset, Item) ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),unsigned short*)
#define ArrayListSetUnsignedIntPtr(ArrayList, ItemOffset, Item)   ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),unsigned int*)
#define ArrayListSetUnsignedLongPtr(ArrayList, ItemOffset, Item)  ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),unsigned long*)
#define ArrayListSetInt64Ptr(ArrayList, ItemOffset, Item)         ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),Int64_t*)
#define ArrayListSetCharPtr(ArrayList, ItemOffset, Item)          ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),char*)
#define ArrayListSetShortPtr(ArrayList, ItemOffset, Item)         ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),short*)
#define ArrayListSetIntPtr(ArrayList, ItemOffset, Item)           ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),int*)
#define ArrayListSetLongPtr(ArrayList, ItemOffset, Item)          ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),long*)
#define ArrayListSetFloatPtr(ArrayList, ItemOffset, Item)         ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),float*)
#define ArrayListSetDoublePtr(ArrayList, ItemOffset, Item)        ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),double*)
#define ArrayListSetLgIndexPtr(ArrayList, ItemOffset, Item)       ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),LgIndex_t*)
#define ArrayListSetEntIndexPtr(ArrayList, ItemOffset, Item)      ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),EntIndex_t*)
#define ArrayListSetSmIntegerPtr(ArrayList, ItemOffset, Item)     ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),SmInteger_t*)
#define ArrayListSetBooleanPtr(ArrayList, ItemOffset, Item)       ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),Boolean_t*)
#define ArrayListSetArbParamPtr(ArrayList, ItemOffset, Item)      ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),ArbParam_t*)
#define ArrayListSetVoidPtr(ArrayList, ItemOffset, Item)          ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),void*)
#define ArrayListSetFunctionPtr(ArrayList, ItemOffset, Item)      ArrayListSetTypedItem((ArrayList),(ItemOffset),(Item),(**)(void))

#define ArrayListAppendUnsignedChar(ArrayList, Item)     ArrayListAppendTypedItem((ArrayList),(Item),unsigned char)
#define ArrayListAppendUnsignedShort(ArrayList, Item)    ArrayListAppendTypedItem((ArrayList),(Item),unsigned short)
#define ArrayListAppendUnsignedInt(ArrayList, Item)      ArrayListAppendTypedItem((ArrayList),(Item),unsigned int)
#define ArrayListAppendUnsignedLong(ArrayList, Item)     ArrayListAppendTypedItem((ArrayList),(Item),unsigned long)
#define ArrayListAppendInt64(ArrayList, Item)            ArrayListAppendTypedItem((ArrayList),(Item),Int64_t)
#define ArrayListAppendChar(ArrayList, Item)             ArrayListAppendTypedItem((ArrayList),(Item),char)
#define ArrayListAppendShort(ArrayList, Item)            ArrayListAppendTypedItem((ArrayList),(Item),short)
#define ArrayListAppendInt(ArrayList, Item)              ArrayListAppendTypedItem((ArrayList),(Item),int)
#define ArrayListAppendLong(ArrayList, Item)             ArrayListAppendTypedItem((ArrayList),(Item),long)
#define ArrayListAppendFloat(ArrayList, Item)            ArrayListAppendTypedItem((ArrayList),(Item),float)
#define ArrayListAppendDouble(ArrayList, Item)           ArrayListAppendTypedItem((ArrayList),(Item),double)
#define ArrayListAppendLgIndex(ArrayList, Item)          ArrayListAppendTypedItem((ArrayList),(Item),LgIndex_t)
#define ArrayListAppendEntIndex(ArrayList, Item)         ArrayListAppendTypedItem((ArrayList),(Item),EntIndex_t)
#define ArrayListAppendSmInteger(ArrayList, Item)        ArrayListAppendTypedItem((ArrayList),(Item),SmInteger_t)
#define ArrayListAppendBoolean(ArrayList, Item)          ArrayListAppendTypedItem((ArrayList),(Item),Boolean_t)
#define ArrayListAppendArbParam(ArrayList, Item)         ArrayListAppendTypedItem((ArrayList),(Item),ArbParam_t)
#define ArrayListAppendUnsignedCharPtr(ArrayList, Item)  ArrayListAppendTypedItem((ArrayList),(Item),unsigned char*)
#define ArrayListAppendUnsignedShortPtr(ArrayList, Item) ArrayListAppendTypedItem((ArrayList),(Item),unsigned short*)
#define ArrayListAppendUnsignedIntPtr(ArrayList, Item)   ArrayListAppendTypedItem((ArrayList),(Item),unsigned int*)
#define ArrayListAppendUnsignedLongPtr(ArrayList, Item)  ArrayListAppendTypedItem((ArrayList),(Item),unsigned long*)
#define ArrayListAppendInt64Ptr(ArrayList, Item)         ArrayListAppendTypedItem((ArrayList),(Item),Int64_t*)
#define ArrayListAppendCharPtr(ArrayList, Item)          ArrayListAppendTypedItem((ArrayList),(Item),char*)
#define ArrayListAppendShortPtr(ArrayList, Item)         ArrayListAppendTypedItem((ArrayList),(Item),short*)
#define ArrayListAppendIntPtr(ArrayList, Item)           ArrayListAppendTypedItem((ArrayList),(Item),int*)
#define ArrayListAppendLongPtr(ArrayList, Item)          ArrayListAppendTypedItem((ArrayList),(Item),long*)
#define ArrayListAppendFloatPtr(ArrayList, Item)         ArrayListAppendTypedItem((ArrayList),(Item),float*)
#define ArrayListAppendDoublePtr(ArrayList, Item)        ArrayListAppendTypedItem((ArrayList),(Item),double*)
#define ArrayListAppendLgIndexPtr(ArrayList, Item)       ArrayListAppendTypedItem((ArrayList),(Item),LgIndex_t*)
#define ArrayListAppendEntIndexPtr(ArrayList, Item)      ArrayListAppendTypedItem((ArrayList),(Item),EntIndex_t*)
#define ArrayListAppendSmIntegerPtr(ArrayList, Item)     ArrayListAppendTypedItem((ArrayList),(Item),SmInteger_t*)
#define ArrayListAppendBooleanPtr(ArrayList, Item)       ArrayListAppendTypedItem((ArrayList),(Item),Boolean_t*)
#define ArrayListAppendArbParamPtr(ArrayList, Item)      ArrayListAppendTypedItem((ArrayList),(Item),ArbParam_t*)
#define ArrayListAppendVoidPtr(ArrayList, Item)          ArrayListAppendTypedItem((ArrayList),(Item),void*)
#define ArrayListAppendFunctionPtr(ArrayList, Item)      ArrayListAppendTypedItem((ArrayList),(Item),(**)(void))

#endif /* ARRLIST_h */
