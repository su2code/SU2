/*
 *****************************************************************
 *****************************************************************
 *******                                                  ********
 ****** Copyright (C) 1988-2010 Tecplot, Inc.              *******
 *******                                                  ********
 *****************************************************************
 *****************************************************************
 */
#if !defined AUXDATA_h
#define AUXDATA_h

#if defined EXTERN
#  undef EXTERN
#endif
#if defined AUXDATAMODULE
#  define EXTERN
#else
#  define EXTERN extern
#endif

/*
 * For building pltview.exe under Windows, we use
 * tecio.dll (which is linked to pltview).
 * Since pltview.exe uses a few of the
 * functions here, they need to be exported into
 * the tecio.dll, thus "TECXXX.h" is included for the
 * LIBFUNCTION & LIBCALL keywords. They are not
 * documented with the other TECXXX() functions,
 * however.
 *
 * If pltview requires other AuxData functions
 * in the future, they can be added to the dll
 * by adding LIBFUNCTION & LIBCALL as in
 * AuxDataDealloc(), etc. below.
 *
 * When building the tecplot kernel, LIBFUNCTION
 * and LIBCALL are nop's.
 *
 */
#include "TECXXX.h"

/**
 */
EXTERN Boolean_t AuxDataIsValidNameChar(char      Char,
                                        Boolean_t IsLeadChar);
/**
 */
EXTERN Boolean_t AuxDataIsValidName(const char *Name);

/**
 */
EXTERN AuxData_pa AuxDataAlloc(void);

/**
 */
LIBFUNCTION void LIBCALL AuxDataDealloc(AuxData_pa *AuxData);

/**
 */
EXTERN Boolean_t AuxDataItemDestructor(void       *ItemRef,
                                       ArbParam_t  ClientData);
/**
 */
EXTERN AuxData_pa AuxDataCopy(AuxData_pa AuxData,
                              Boolean_t  ConsiderRetain);

/**
 */
LIBFUNCTION LgIndex_t LIBCALL AuxDataGetNumItems(AuxData_pa AuxData);

/**
 */
EXTERN Boolean_t AuxDataGetItemIndex(AuxData_pa AuxData,
                                     const char *Name,
                                     LgIndex_t  *ItemIndex);
/**
 */
LIBFUNCTION void LIBCALL AuxDataGetItemByIndex(AuxData_pa    AuxData,
			                                   LgIndex_t     Index,
			                                   const char    **Name,
			                                   ArbParam_t    *Value,
			                                   AuxDataType_e *Type,
			                                   Boolean_t     *Retain);

/**
 */
EXTERN Boolean_t AuxDataGetItemByName(AuxData_pa    AuxData,
                                      const char    *Name,
                                      ArbParam_t    *Value,
                                      AuxDataType_e *Type,
                                      Boolean_t     *Retain);

/**
 */
EXTERN Boolean_t AuxDataGetBooleanItemByName(AuxData_pa     AuxData,
                                             const char    *Name,
                                             Boolean_t     *Value,
                                             AuxDataType_e *Type,
                                             Boolean_t     *Retain);

/**
 */
EXTERN Boolean_t AuxDataSetItem(AuxData_pa    AuxData,
                                const char    *Name,
                                ArbParam_t    Value,
                                AuxDataType_e Type,
                                Boolean_t     Retain);

/**
 */
EXTERN Boolean_t AuxDataDeleteItemByName(AuxData_pa AuxData,
                                         const char *Name);

/**
 */
EXTERN Boolean_t AuxDataAppendItems(AuxData_pa TargetAuxData,
                                    AuxData_pa SourceAuxData);
/**
 */
EXTERN void AuxDataDeleteItemByIndex(AuxData_pa AuxData,
                                     LgIndex_t  Index);

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

#endif /* !defined AUXDATA_h */
