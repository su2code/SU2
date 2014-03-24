#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE


/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/



#define SETMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "ALLOC.h"
#include "SET.h"

/* * SET FUNCTIONS * */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if InitNumZones > InitNumVars
#else
#endif
#if ZoneExpansionFactor > VarExpansionFactor
#else
#endif
#else
#define SetInitSize          (PadOut(1,SetBitSize))
#define SetExpansionFactor   2
#endif

using tecplot::strutil::translate;

/*
 */
Set_pa AllocSet(Boolean_t show_error_msg)
{
    Set_pa Set = ALLOC_ITEM(struct _Set_a, "Set header");
    if (Set)
    {
        Set->size = SetInitSize;
        Set->data = ALLOC_ARRAY(SetInitSize / SetBitSize, SetData_t, "Set data");
        if (Set->data == NULL)
            DeallocSet(&Set);
        else
            ClearSet(Set);
    }
    if ((Set == NULL) && show_error_msg)
    {
#     if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#     else
        fprintf(stderr, "Out of memory for sets");
#     endif
    }
    return Set;
} /* AllocSet() */


/*
 */
void DeallocSet(Set_pa *Set)
{
    if (Set && *Set)
    {
        if ((*Set)->data)
            FREE_ARRAY((*Set)->data, "Set data");
        FREE_ITEM(*Set, "Set header");
        *Set = NULL;
    }
} /* DeallocSet() */


/**
 * This function adapts the DeallocSet function to work with the
 * ArrayList's deallocation callback.
 */
Boolean_t SetItemDestructor(void       *ItemRef,
                            ArbParam_t  ClientData)
{
    Set_pa *SetRef = (Set_pa *)ItemRef;

    REQUIRE(VALID_REF(SetRef));
    REQUIRE(VALID_REF(*SetRef) || *SetRef == NULL);
    UNUSED(ClientData);

    if (*SetRef != NULL)
        DeallocSet(SetRef);

    ENSURE(*SetRef == NULL);
    return TRUE;
}


/*
 */
Boolean_t ExpandSet(Set_pa     Set,
                    SetIndex_t max_val,
                    Boolean_t  show_error_msg)
{
    SetData_t  *data;
    SetIndex_t  new_size;

    REQUIRE(max_val >= 0);

    if (!Set)
    {
        if (show_error_msg)
        {
#         if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#         else
            fprintf(stderr, "Null Set expand");
#         endif
        }
        return FALSE;
    }

    if (max_val <= Set->size)
        return TRUE;

    new_size = Set->size;
    while (new_size < max_val)
        new_size *= SetExpansionFactor;

    new_size = PadOut(new_size, SetBitSize);

    data = ALLOC_ARRAY(new_size / SetBitSize, SetData_t, "new Set data");

    if (!data)
    {
        if (show_error_msg)
        {
#         if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#         else
            fprintf(stderr, "Out of memory for sets");
#         endif
        }
        return FALSE;
    }
    size_t old_set_size_in_bytes = sizeof(data[0]) * (Set->size / SetBitSize);
    memcpy(data, Set->data, old_set_size_in_bytes);

    size_t new_set_size_in_bytes = sizeof(data[0]) * (new_size / SetBitSize);
    size_t numBytesToReset = new_set_size_in_bytes - old_set_size_in_bytes;
    memset(((char*)data) + old_set_size_in_bytes, 0, numBytesToReset);

    FREE_ARRAY(Set->data, "old Set data");
    Set->data = data;
    Set->size = new_size;
    return TRUE;
} /* ExpandSet() */


/*
 */
Boolean_t CopySet(Set_pa    dst,
                  Set_pa    src,
                  Boolean_t show_error_msg)
{
    if (dst && dst->data &&
        src && src->data &&
        ExpandSet(dst, src->size, show_error_msg))
    {
        SetIndex_t src_size_in_words = src->size / SetBitSize;
        size_t numBytesToCopy = sizeof(dst->data[0]) * src_size_in_words;
        memcpy(dst->data, src->data, numBytesToCopy);

        SetIndex_t dst_size_in_words = dst->size / SetBitSize;
        CHECK(dst_size_in_words>=src_size_in_words); // ...guaranteed by above ExpandSet() call
        size_t numBytesToReset = sizeof(dst->data[0]) * (dst_size_in_words - src_size_in_words);
        memset((char*)(dst->data + src_size_in_words), 0, numBytesToReset);

        return TRUE;
    }
    else
        return FALSE;
} /* CopySet() */


/*
 */
Boolean_t AppendSet(Set_pa dst,
                    Set_pa src,
                    Boolean_t show_error_msg)
{
    if (dst && dst->data &&
        src && src->data)
    {
        SetIndex_t member;
        ForAllMembersInSet(member, src)
        {
            if (!AddToSet(dst, member, show_error_msg))
                return FALSE;
        }
        return TRUE;
    }
    else
        return FALSE;
} /* AppendSet() */


/*
 */
void ClearSet(Set_pa Set)
{
    if (Set && Set->data)
        memset(Set->data, 0, Set->size / SetBitSize * sizeof(Set->data[0]));
} /* ClearSet() */


#if defined USE_FUNCTIONS_FOR_SETS
/*
 */
Boolean_t AddToSet(Set_pa     Set,
                   SetIndex_t member,
                   Boolean_t  show_error_msg)
{
    REQUIRE(member >= 0);
    if (Set &&
        Set->data &&
        ((member + 1 <= Set->size) || ExpandSet(Set, member + 1, show_error_msg)))
    {
        SetIndex_t word = member / SetBitSize;
        SetData_t  bit = (SetData_t)1 << (member % SetBitSize);
        Set->data[word] |= bit;
        return TRUE;
    }
    else
        return FALSE;
} /* AddToSet() */
#endif


/*
 */
void RemoveFromSet(Set_pa     Set,
                   SetIndex_t member)
{
    REQUIRE(member >= 0);
    if (Set && (member < Set->size) && Set->data)
    {
        SetIndex_t word = member / SetBitSize;
        SetData_t  bit = (SetData_t)1 << (member % SetBitSize);
        Set->data[word] &= (((SetData_t) - 1) ^ bit);
    }
} /* RemoveFromSet() */


/**
 * Similar to RemoveFromSet except it shifts the Set.
 */
void DeleteSetMember(Set_pa     Set,
                     SetIndex_t Member)
{
    SetIndex_t LastMember;

    REQUIRE(VALID_REF(Set));
    REQUIRE(Member >= 0);

    LastMember = GetPrevMember(Set, BAD_SET_VALUE);
    if (Member <= LastMember)
    {
        ShiftSet(Set, Member + 1, LastMember, -1);
        RemoveFromSet(Set, LastMember);
    }
}


/**
 * Similar to AddToSet except that if the new member is within the currently
 * defined set the members are shifted accordingly.
 */
Boolean_t InsertSetMember(Set_pa     Set,
                          SetIndex_t Member,
                          Boolean_t  ShowErrMsg)
{
    Boolean_t  IsOk = TRUE;
    SetIndex_t OrigLastMember;

    REQUIRE(VALID_REF(Set));

    /* first, determine if we need to shift the set */
    OrigLastMember = GetPrevMember(Set, BAD_SET_VALUE);
    if (Member <= OrigLastMember)
    {
        IsOk = ExpandSet(Set, (OrigLastMember + 1) + 1, ShowErrMsg);
        ShiftSet(Set, Member, OrigLastMember, 1);
    }

    if (IsOk)
        IsOk = AddToSet(Set, Member, ShowErrMsg);

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}

#if defined USE_FUNCTIONS_FOR_SETS
/*
 */
Boolean_t InSet(Set_pa     Set,
                SetIndex_t member)
{
    /*
     * Sometimes InSet is called with negative numbers.  This is not correct, but
     * its what we have to work with.  Maybe some day, we can make this assertion.
    REQUIRE(member>=0);
     */
    if (Set && (0 <= member && member < Set->size))
    {
        SetIndex_t word = member / SetBitSize;
        SetData_t  bit = (SetData_t)1 << (member % SetBitSize);
        return (Set->data[word]&bit) != 0;
    }
    else
        return FALSE;
} /* InSet() */
#endif


/*
 */
Boolean_t IsEmpty(Set_pa Set)
{
    if (Set && Set->data)
    {
        SetIndex_t set_size_in_words = Set->size / SetBitSize;
        SetIndex_t word;
        for (word = 0; word < set_size_in_words; word++)
            if (Set->data[word] != 0)
                return FALSE;
    }
    return TRUE;
} /* IsEmpty() */


/*
 */
Boolean_t HasVoids(Set_pa Set)
{
    Boolean_t  Result = FALSE;
    SetIndex_t ContiguousMember = 0;
    SetIndex_t Member = 0;

    REQUIRE(VALID_REF(Set));

    /* look for voids in the set */
    ForAllMembersInSet(Member, Set)
    {
        if (Member == ContiguousMember)
        {
            ContiguousMember++;
        }
        else
        {
            Result = TRUE;
            break;
        }
    }

    ENSURE(VALID_BOOLEAN(Result));
    return Result;
}


/*
 */
SetIndex_t MemberCount(Set_pa Set)
{
    SetIndex_t count = 0;
    if (Set && Set->data)
    {
        SetIndex_t set_size_in_words = Set->size / SetBitSize;
        SetIndex_t word;
        for (word = 0; word < set_size_in_words; word++)
        {
            SetData_t word_val = Set->data[word];
            while (word_val)
            {
                if (word_val&1)
                    count++;
                word_val = word_val >> 1;
            }
        }
    }
    return count;
} /* MemberCount() */


/*
 */
SetIndex_t GetNextMember(Set_pa     Set,
                         SetIndex_t start_at)
{
    SetIndex_t next_member = BAD_SET_VALUE;
    if (Set && Set->data)
    {
        SetIndex_t set_size_in_words = Set->size / SetBitSize;
        SetIndex_t word;
        SetData_t  word_val = 0;
        int bit;
        if (start_at == BAD_SET_VALUE)
        {
            word = 0;
            bit = 0;
            if (word < set_size_in_words)
                word_val = Set->data[0];
        }
        else if (start_at + 1 < Set->size)
        {
            word = (start_at + 1) / SetBitSize;
            bit = static_cast<int>((start_at + 1) % SetBitSize);
            if (word < set_size_in_words)
                word_val = Set->data[word] >> bit;
        }
        else
        {
            return BAD_SET_VALUE;
        }
        while ((word < set_size_in_words) && (word_val == 0))
        {
            word++;
            bit = 0;
            if (word < set_size_in_words)
                word_val = Set->data[word];
        }
        if (word < set_size_in_words)
        {
            while (!(word_val&1))
            {
                word_val >>= 1;
                bit++;
            }
            next_member = word * SetBitSize + bit;
        }
    }
    return next_member;
} /* GetNextMember() */


/*
 */
SetIndex_t GetPrevMember(Set_pa     Set,
                         SetIndex_t start_at)
{
    SetIndex_t next_member = BAD_SET_VALUE;
    if (Set && Set->data)
    {
        SetIndex_t set_size_in_words = Set->size / SetBitSize;
        SetIndex_t word;
        SetData_t  word_val = 0;
        int bit;
        if (start_at == BAD_SET_VALUE)
        {
            word = set_size_in_words - 1;
            bit = SetBitSize - 1;
            if (word >= 0)
                word_val = Set->data[word];
        }
        else if (start_at > 0)
        {
            word = (start_at - 1) / SetBitSize;
            bit = static_cast<int>((start_at - 1) % SetBitSize);
            if (word >= 0)
                word_val = Set->data[word] << (SetBitSize - bit - 1);
        }
        else
        {
            return BAD_SET_VALUE;
        }
        while ((word >= 0) && (word_val == 0))
        {
            word--;
            bit = static_cast<int>(SetBitSize - 1);
            if (word >= 0)
                word_val = Set->data[word] << (SetBitSize - bit - 1);
        }
        if (word >= 0)
        {
            while (!(word_val&SetLastBit))
            {
                word_val <<= 1;
                bit--;
            }
            next_member = word * SetBitSize + bit;
        }
    }
    return next_member;
} /* GetPrevMember() */


/*
 */
Boolean_t EqualSets(Set_pa set1,
                    Set_pa set2)
{
    SetIndex_t set1_size_in_words,
    set2_size_in_words,
    min_set_size_in_words,
    ii;
    if (!set1 || !set2)
        return FALSE;

    set1_size_in_words = set1->size / SetBitSize;
    set2_size_in_words = set2->size / SetBitSize;
    min_set_size_in_words = MIN(set1_size_in_words, set2_size_in_words);

    for (ii = 0; ii < min_set_size_in_words; ii++)
        if (set1->data[ii] != set2->data[ii])
            return FALSE;
    for (ii = min_set_size_in_words; ii < set1_size_in_words; ii++)
        if (set1->data[ii] != 0)
            return FALSE;
    for (ii = min_set_size_in_words; ii < set2_size_in_words; ii++)
        if (set2->data[ii] != 0)
            return FALSE;

    return TRUE;

} /* EqualSets() */


Boolean_t IsSubSet(Set_pa childset,
                   Set_pa parentset)
{
    SetIndex_t s;

    ForAllMembersInSet(s, childset)
    {
        if (!InSet(parentset, s))
            return (FALSE);
    }

    return (TRUE);

} /* IsSubSet() */






/*
 *   functions added 11/7 by byron.  These are roughed in for now and could
 *   stand to be optimized later.....
 */

/*
 *  Return the number of members in a set that preceed a given member.
 */

SetIndex_t MemberOffset(Set_pa     Set,
                        SetIndex_t Member)
{
    SetIndex_t I;
    SetIndex_t Offset = -1;
    if (InSet(Set, Member))
    {
        for (I = 0; I <= Member; I++)
        {
            if (InSet(Set, I))
                Offset++;
        }
    }
    return (Offset);
}

/*
 *  Return the position in the set of the nth member of a set.
 */
SetIndex_t OffsetMember(Set_pa     Set,
                        SetIndex_t Offset)
{
    SetIndex_t I;
    SetIndex_t Member = BAD_SET_VALUE;
    for (I = 0; I <= Offset; I++)
    {
        Member = GetNextMember(Set, Member);
        if (Member == BAD_SET_VALUE)
            break;
    }
    return (Member);
}

Boolean_t CopySetMember(Set_pa     DstSet,
                        SetIndex_t DstOffset,
                        Set_pa     SrcSet,
                        SetIndex_t SrcOffset)
{
    if (InSet(SrcSet, SrcOffset))
        return (AddToSet(DstSet, DstOffset, TRUE));
    else
        RemoveFromSet(DstSet, DstOffset);
    return (TRUE);
}



/*
 *  Initial:
 *                v---ShiftPos1   v--ShiftPos2
 * +-------------------------------------+
 * | | | | | | | |x| | | | | | | |x| | | |
 * +-------------------------------------+
 *
 * Shift +2
 *                    v---ShiftPos1   v--ShiftPos2
 * +-------------------------------------+
 * | | | | | | | | | |x| | | | | | | |x| |
 * +-------------------------------------+
 *
 *
 * Shift all bits between ShiftPos1 and ShiftPos2
 * by ShiftAmount.  The bits that the shift
 * replaces fill in the hole left by the shift
 *
 *
 */
void ShiftSet(Set_pa     Set,
              SetIndex_t ShiftPos1,
              SetIndex_t ShiftPos2,
              SetIndex_t ShiftAmount)
{
    Set_pa     NewSet;
    SetIndex_t DPos;
    SetIndex_t SPos;

    if ((Set == NULL) || (IsEmpty(Set)))
        return;

    NewSet = AllocSet(TRUE);

    if (NewSet == NULL)
        return;

    if (!CopySet(NewSet, Set, TRUE))
        return;

    if (ShiftAmount < 0)
    {
        DPos = ShiftPos2;
        SPos = ShiftPos1 - 1;
        while (DPos > ShiftPos2 + ShiftAmount)
            CopySetMember(NewSet, DPos--, Set, SPos--);
        SPos = ShiftPos2;
        while (SPos >= ShiftPos1)
            CopySetMember(NewSet, DPos--, Set, SPos--);
    }
    else if (ShiftAmount > 0)
    {
        DPos       = ShiftPos1;
        SPos       = ShiftPos2 + 1;
        while (DPos < ShiftPos1 + ShiftAmount)
            CopySetMember(NewSet, DPos++, Set, SPos++);
        SPos = ShiftPos1;
        while (SPos <= ShiftPos2)
            CopySetMember(NewSet, DPos++, Set, SPos++);
    }
    CopySet(Set, NewSet, TRUE);
    DeallocSet(&NewSet);
}
