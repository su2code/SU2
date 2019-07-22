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
typedef struct ___135* ___134;
 #endif
typedef enum { ArrayListType_UInt8, ArrayListType_UInt16, ArrayListType_UInt32, ArrayListType_UInt64, ArrayListType_Int64, ArrayListType_Char, ArrayListType_Short, ArrayListType_Int, ArrayListType_Long, ArrayListType_Float, ArrayListType_Double, ArrayListType_LgIndex, ArrayListType_EntIndex, ArrayListType_SmInteger, ArrayListType_Boolean, ArrayListType_ArbParam, ArrayListType_UInt8Ptr, ArrayListType_UInt16Ptr, ArrayListType_UInt32Ptr, ArrayListType_UInt64Ptr, ArrayListType_Int64Ptr, ArrayListType_CharPtr, ArrayListType_ShortPtr, ArrayListType_IntPtr, ArrayListType_LongPtr, ArrayListType_FloatPtr, ArrayListType_DoublePtr, ArrayListType_LgIndexPtr, ArrayListType_EntIndexPtr, ArrayListType_SmIntegerPtr, ArrayListType_BooleanPtr, ArrayListType_ArbParamPtr, ArrayListType_VoidPtr, ArrayListType_FunctionPtr, ArrayListType_Any, END_ArrayListType_e, ArrayListType_Invalid = ___329 } ArrayListType_e; typedef union { uint8_t         UInt8; uint16_t        UInt16; uint32_t        UInt32; uint64_t        UInt64; int64_t         ___1969; char            ___472; short           Short; int             Int; long            Long; float           Float; double          Double; ___2227       ___2225; ___1172      ___1170; int32_t     ___3634; ___372       ___350; ___90      ___88; uint8_t*        UInt8Ptr; uint16_t*       UInt16Ptr; uint32_t*       UInt32Ptr; uint64_t*       UInt64Ptr; int64_t*        ___1970; char*           ___474; short*          ___3561; int*            ___1986; long*           ___2321; float*          ___1436; double*         ___1110; ___2227*      ___2226; ___1172*     ___1171; int32_t*    ___3635; ___372*      ___371; ___90*     ___89; void*           ___4440; void (*___1541)(void); } ArrayListItem_u; typedef ___372 (*ArrayListItemVisitor_pf)(void*      ___2098, ___90 ___494);
 #if 0 
{ REQUIRE(VALID_REF(___4239)); REQUIRE(VALID_REF(*___4239) || *___4239 == NULL); ___372 ___1096 = ___4226; <type>* ___4239 = static_cast<<type>*>(___2098); ENSURE(VALID_BOOLEAN(___1096)); return ___1096; }
 #endif
typedef ArrayListItemVisitor_pf ArrayListItemDestructor_pf; typedef ___372 (*ArrayListItemDuplicator_pf)(void*      ___3949, void*      ___3645, ___90 ___494);
 #if 0 
{ REQUIRE(VALID_REF(___3950)); REQUIRE(VALID_REF(___3646)); REQUIRE(VALID_REF(*___3646) || *___3646 == NULL); ___372 ___2040 = ___4226; <type>* ___3950 = static_cast<<type>*>(___3949); <type>* ___3646 = static_cast<<type>*>(___3645); ENSURE(VALID_BOOLEAN(___2040)); return ___2040; }
 #endif
typedef ___2227 (*ArrayListCapacityRequestAdjuster_pf)(___134 ___94, ___2227    ___693, ___2227    ___3354, ___90   ___494);
 #if 0 
{ REQUIRE(ArrayListIsValid(___94)); REQUIRE((___3354 == 0 && ___693 == 0) || ___3354 > ___94->___439); ___2227 ___3359; ENSURE(___3359 == 0 || ___3359 >= ___3354); return ___3359; }
 #endif
struct ___135 { char*            Array; ArrayListType_e  ___4236; int32_t      ___2102; ___2227        ___684; ___2227        ___439; ___372        ___2080; ArrayListCapacityRequestAdjuster_pf ___440; ___90                          ___441; }; typedef int (STDCALL *ArrayListItemComparator_pf)(ArrayListItem_u ___2087, ArrayListItem_u ___2088, ___90      ___494); EXTERN ___372 ArrayListIsValid(___134 ___94); EXTERN ArrayListType_e ArrayListGetType(___134 ___94); EXTERN ___372 ArrayListEnlargeCapacity(___134 ___94, ___2227    ___3354); EXTERN ___134 ArrayListAlloc(___2227                           ___1188, ArrayListType_e                     ___4236, ArrayListCapacityRequestAdjuster_pf ___440 = 0, ___90                          ___441 = 0); EXTERN void ArrayListDealloc(___134*              ___94, ArrayListItemDestructor_pf ___2094 = 0, ___90                 ___494 = 0); EXTERN void ArrayListDeleteAllItems(___134               ___94, ArrayListItemDestructor_pf ___2094 = 0, ___90                 ___494 = 0); EXTERN void ArrayListDeleteItems(___134               ___94, ___2227                  ___2097, ___2227                  ___684, ArrayListItemDestructor_pf ___2094 = 0, ___90                 ___494 = 0); EXTERN void ArrayListDeleteItem(___134               ___94, ___2227                  ___2097, ArrayListItemDestructor_pf ___2094 = 0, ___90                 ___494 = 0); EXTERN ___134 ArrayListRemoveItems(___134 ___94, ___2227    ___2097, ___2227    ___684); EXTERN ArrayListItem_u ArrayListRemoveItem(___134 ___94, ___2227    ___2097); EXTERN ___372 ArrayListInsertItem(___134    ___94, ___2227       ___2097, ArrayListItem_u ___2086); EXTERN ___372 ArrayListInsert(___134 ___3946, ___2227    ___2097, ___134 ___3642); EXTERN ___372 ArrayListVisitItems(___134            ___94, ___2227               ___2097, ___2227               ___684, ArrayListItemVisitor_pf ___2103, ___90              ___494); EXTERN ___134 ArrayListGetItems(___134 ___94, ___2227    ___2097, ___2227    ___684); EXTERN ArrayListItem_u ArrayListGetItem(___134 ___94, ___2227    ___2097); EXTERN ___372 ArrayListSetItem(___134               ___94, ___2227                  ___2097, ArrayListItem_u            ___2086, ArrayListItemDestructor_pf ___2094 = 0, ___90                 ___494 = 0); EXTERN ___372 ArrayListAppendItem(___134    ___94, ArrayListItem_u ___2086); EXTERN ___372 ArrayListAppend(___134 ___3946, ___134 ___3642); EXTERN ___134 ArrayListCopy(___134               ___94, ArrayListItemDuplicator_pf ___2095 = 0, ___90                 ___494 = 0); EXTERN void* ArrayListToArray(___134               ___94, ArrayListItemDuplicator_pf ___2095, ___90                 ___494); EXTERN ___134 ArrayListFromArray(void*                      ___3642, ___2227                  ___684, ArrayListType_e            ___4236, ArrayListItemDuplicator_pf ___2095 = 0, ___90                 ___494 = 0); EXTERN ___372 ArrayListBSearch(___134               ___94, ArrayListItem_u            ___2086, ArrayListItemComparator_pf ___535, ___90                 ___494, ___2227*                 ___2096 = 0);
 #if defined USE_MACROS_FOR_FUNCTIONS
 #  define ___112     ___113
 #  define ___115 ArrayListGetItemInternalRef_MACRO
 #  define ___101           ArrayListGetCount_MACRO
 #  define ___131(___94, ___2097)            ___124(___94, ___2097, uint8_t)
 #  define ___125(___94, ___2097)           ___124(___94, ___2097, uint16_t)
 #  define ___127(___94, ___2097)           ___124(___94, ___2097, uint32_t)
 #  define ___129(___94, ___2097)           ___124(___94, ___2097, uint64_t)
 #  define ___110(___94, ___2097)            ___124(___94, ___2097, int64_t)
 #  define ___99(___94, ___2097)             ___124(___94, ___2097, char)
 #  define ___120(___94, ___2097)            ___124(___94, ___2097, short)
 #  define ___109(___94, ___2097)              ___124(___94, ___2097, int)
 #  define ___118(___94, ___2097)             ___124(___94, ___2097, long)
 #  define ___106(___94, ___2097)            ___124(___94, ___2097, float)
 #  define ___102(___94, ___2097)           ___124(___94, ___2097, double)
 #  define ___116(___94, ___2097)          ___124(___94, ___2097, ___2227)
 #  define ___104(___94, ___2097)         ___124(___94, ___2097, ___1172)
 #  define ___122(___94, ___2097)        ___124(___94, ___2097, int32_t)
 #  define ___97(___94, ___2097)          ___124(___94, ___2097, ___372)
 #  define ___95(___94, ___2097)         ___124(___94, ___2097, ___90)
 #  define ___132(___94, ___2097)         ___124(___94, ___2097, uint8_t*)
 #  define ___126(___94, ___2097)        ___124(___94, ___2097, uint16_t*)
 #  define ___128(___94, ___2097)        ___124(___94, ___2097, uint32_t*)
 #  define ___130(___94, ___2097)        ___124(___94, ___2097, uint64_t*)
 #  define ___111(___94, ___2097)         ___124(___94, ___2097, int64_t*)
 #  define ___100(___94, ___2097)          ___124(___94, ___2097, char*)
 #  define ___121(___94, ___2097)         ___124(___94, ___2097, short*)
 #  define ___114(___94, ___2097)           ___124(___94, ___2097, int*)
 #  define ___119(___94, ___2097)          ___124(___94, ___2097, long*)
 #  define ___107(___94, ___2097)         ___124(___94, ___2097, float*)
 #  define ___103(___94, ___2097)        ___124(___94, ___2097, double*)
 #  define ___117(___94, ___2097)       ___124(___94, ___2097, ___2227*)
 #  define ___105(___94, ___2097)      ___124(___94, ___2097, ___1172*)
 #  define ___123(___94, ___2097)     ___124(___94, ___2097, int32_t*)
 #  define ___98(___94, ___2097)       ___124(___94, ___2097, ___372*)
 #  define ___96(___94, ___2097)      ___124(___94, ___2097, ___90*)
 #  define ___133(___94, ___2097)          ___124(___94, ___2097, void*)
 #  define ___108(___94, ___2097)      ___124(___94, ___2097, (**)(void))
 #else
 #  define ___112     ArrayListGetInternalRef_FUNC
 #  define ___115 ArrayListGetItemInternalRef_FUNC
 #  define ___101           ArrayListGetCount_FUNC
 #  define ___131(___94, ___2097)                    (*(static_cast<uint8_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___125(___94, ___2097)                  (*(static_cast<uint16_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___127(___94, ___2097)                  (*(static_cast<uint32_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___129(___94, ___2097)                  (*(static_cast<uint64_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___110(___94, ___2097)                    (*(static_cast<int64_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___99(___94, ___2097)                        (*(static_cast<char*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___120(___94, ___2097)                      (*(static_cast<short*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___109(___94, ___2097)                          (*(static_cast<int*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___118(___94, ___2097)                        (*(static_cast<long*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___106(___94, ___2097)                      (*(static_cast<float*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___102(___94, ___2097)                    (*(static_cast<double*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___116(___94, ___2097)                (*(static_cast<___2227*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___104(___94, ___2097)              (*(static_cast<___1172*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___122(___94, ___2097)            (*(static_cast<int32_t*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___97(___94, ___2097)                (*(static_cast<___372*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___95(___94, ___2097)              (*(static_cast<___90*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___132(___94, ___2097)                (*(static_cast<uint8_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ArrayListGetUInt16tPtr(___94, ___2097)             (*(static_cast<uint16_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ArrayListGetUInt32tr(___94, ___2097)               (*(static_cast<uint32_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___130(___94, ___2097)              (*(static_cast<uint64_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___111(___94, ___2097)                (*(static_cast<int64_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___100(___94, ___2097)                    (*(static_cast<char**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___121(___94, ___2097)                  (*(static_cast<short**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___114(___94, ___2097)                      (*(static_cast<int**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___119(___94, ___2097)                    (*(static_cast<long**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___107(___94, ___2097)                  (*(static_cast<float**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___103(___94, ___2097)                (*(static_cast<double**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___117(___94, ___2097)            (*(static_cast<___2227**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___105(___94, ___2097)          (*(static_cast<___1172**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___123(___94, ___2097)        (*(static_cast<int32_t**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___98(___94, ___2097)            (*(static_cast<___372**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___96(___94, ___2097)          (*(static_cast<___90**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___133(___94, ___2097)                    (*(static_cast<void**>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #  define ___108(___94, ___2097)             (*(static_cast<**(void)*>(const_cast<void*>(ArrayListGetItemInternalRef_FUNC(___94, ___2097)))))
 #endif
 #if !defined USE_MACROS_FOR_FUNCTIONS
EXTERN void const* ArrayListGetInternalRef_FUNC(___134 ___94); EXTERN void const* ArrayListGetItemInternalRef_FUNC(___134 ___94, ___2227    ___2097); EXTERN ___2227 ArrayListGetCount_FUNC(___134 ___94);
 #endif
 #define ___113(___94)                 static_cast<void const*>((___94)->Array)
 #define ArrayListGetItemInternalRef_MACRO(___94, ___2097) static_cast<void const*>(&((___94)->Array[(___2097)*(___94)->___2102]))
 #define ArrayListGetCount_MACRO(___94)                       ((___94)->___684)
 #define ArrayListGetTypedArrayRef(___94, ___2688)         reinterpret_cast<___2688*>((___94)->Array)
 #define ___124(___94, ___2097, ___2688) (ArrayListGetTypedArrayRef(___94,___2688)[___2097])
 #if defined NO_ASSERTS
 # define ArrayListOffsetWithinCapacity(___94, ___2097) ((___2097) < (___94)->___439)
 #else
 # define ArrayListOffsetWithinCapacity(___94, ___2097) ((assert((___2097) >= 0),___4226) && ((___2097) < (___94)->___439))
 #endif
 #define ___161(___94, ___2097, ___2086, ___2688) \
 ((ArrayListOffsetWithinCapacity((___94), (___2097)) || \
 ArrayListEnlargeCapacity((___94), (___2097)+1)) \
 ? ((void)((ArrayListGetTypedArrayRef((___94),___2688)[(___2097)]) = (___2086)), \
 (((___2097)+1 > (___94)->___684) \
 ? (((___94)->___684 = (___2097)+1), ___4226) \
 : (___4226))) \
 : (___1305))
 #define ArrayListAppendTypedItem(___94, ___2086, ___2688) \
 ((ArrayListOffsetWithinCapacity((___94), (___94)->___684) || \
 ArrayListEnlargeCapacity((___94), (___94)->___684+1)) \
 ? ((void)((ArrayListGetTypedArrayRef((___94),___2688)[(___94)->___684]) = (___2086)), \
 (((___94)->___684 = (___94)->___684+1), ___4226)) \
 : (___1305))
 #define ___168(___94, ___2097, ___2086)            ___161((___94),(___2097),(___2086),uint8_t)
 #define ___162(___94, ___2097, ___2086)           ___161((___94),(___2097),(___2086),uint16_t)
 #define ___164(___94, ___2097, ___2086)           ___161((___94),(___2097),(___2086),uint32_t)
 #define ___166(___94, ___2097, ___2086)           ___161((___94),(___2097),(___2086),uint64_t)
 #define ___150(___94, ___2097, ___2086)            ___161((___94),(___2097),(___2086),int64_t)
 #define ___140(___94, ___2097, ___2086)             ___161((___94),(___2097),(___2086),char)
 #define ___157(___94, ___2097, ___2086)            ___161((___94),(___2097),(___2086),short)
 #define ___149(___94, ___2097, ___2086)              ___161((___94),(___2097),(___2086),int)
 #define ___155(___94, ___2097, ___2086)             ___161((___94),(___2097),(___2086),long)
 #define ___146(___94, ___2097, ___2086)            ___161((___94),(___2097),(___2086),float)
 #define ___142(___94, ___2097, ___2086)           ___161((___94),(___2097),(___2086),double)
 #define ___153(___94, ___2097, ___2086)          ___161((___94),(___2097),(___2086),___2227)
 #define ___144(___94, ___2097, ___2086)         ___161((___94),(___2097),(___2086),___1172)
 #define ___159(___94, ___2097, ___2086)        ___161((___94),(___2097),(___2086),int32_t)
 #define ___138(___94, ___2097, ___2086)          ___161((___94),(___2097),(___2086),___372)
 #define ___136(___94, ___2097, ___2086)         ___161((___94),(___2097),(___2086),___90)
 #define ___169(___94, ___2097, ___2086)         ___161((___94),(___2097),(___2086),uint8_t*)
 #define ___163(___94, ___2097, ___2086)       ___161((___94),(___2097),(___2086),uint16_t*)
 #define ___165(___94, ___2097, ___2086)         ___161((___94),(___2097),(___2086),uint32_t*)
 #define ___167(___94, ___2097, ___2086)        ___161((___94),(___2097),(___2086),uint64_t*)
 #define ___151(___94, ___2097, ___2086)         ___161((___94),(___2097),(___2086),int64_t*)
 #define ___141(___94, ___2097, ___2086)          ___161((___94),(___2097),(___2086),char*)
 #define ___158(___94, ___2097, ___2086)         ___161((___94),(___2097),(___2086),short*)
 #define ___152(___94, ___2097, ___2086)           ___161((___94),(___2097),(___2086),int*)
 #define ___156(___94, ___2097, ___2086)          ___161((___94),(___2097),(___2086),long*)
 #define ___147(___94, ___2097, ___2086)         ___161((___94),(___2097),(___2086),float*)
 #define ___143(___94, ___2097, ___2086)        ___161((___94),(___2097),(___2086),double*)
 #define ___154(___94, ___2097, ___2086)       ___161((___94),(___2097),(___2086),___2227*)
 #define ___145(___94, ___2097, ___2086)      ___161((___94),(___2097),(___2086),___1172*)
 #define ___160(___94, ___2097, ___2086)     ___161((___94),(___2097),(___2086),int32_t*)
 #define ___139(___94, ___2097, ___2086)       ___161((___94),(___2097),(___2086),___372*)
 #define ___137(___94, ___2097, ___2086)      ___161((___94),(___2097),(___2086),___90*)
 #define ___170(___94, ___2097, ___2086)          ___161((___94),(___2097),(___2086),void*)
 #define ___148(___94, ___2097, ___2086)      ___161((___94),(___2097),(___2086),(**)(void))
 #define ArrayListAppendUInt8(___94, ___2086)            ArrayListAppendTypedItem((___94),(___2086),uint8_t)
 #define ArrayListAppendUInt16(___94, ___2086)           ArrayListAppendTypedItem((___94),(___2086),uint16_t)
 #define ArrayListAppendUInt32(___94, ___2086)           ArrayListAppendTypedItem((___94),(___2086),uint32_t)
 #define ArrayListAppendUInt64(___94, ___2086)           ArrayListAppendTypedItem((___94),(___2086),uint64_t)
 #define ArrayListAppendInt64(___94, ___2086)            ArrayListAppendTypedItem((___94),(___2086),int64_t)
 #define ArrayListAppendChar(___94, ___2086)             ArrayListAppendTypedItem((___94),(___2086),char)
 #define ArrayListAppendShort(___94, ___2086)            ArrayListAppendTypedItem((___94),(___2086),short)
 #define ArrayListAppendInt(___94, ___2086)              ArrayListAppendTypedItem((___94),(___2086),int)
 #define ArrayListAppendLong(___94, ___2086)             ArrayListAppendTypedItem((___94),(___2086),long)
 #define ArrayListAppendFloat(___94, ___2086)            ArrayListAppendTypedItem((___94),(___2086),float)
 #define ArrayListAppendDouble(___94, ___2086)           ArrayListAppendTypedItem((___94),(___2086),double)
 #define ArrayListAppendLgIndex(___94, ___2086)          ArrayListAppendTypedItem((___94),(___2086),___2227)
 #define ArrayListAppendEntIndex(___94, ___2086)         ArrayListAppendTypedItem((___94),(___2086),___1172)
 #define ArrayListAppendSmInteger(___94, ___2086)        ArrayListAppendTypedItem((___94),(___2086),int32_t)
 #define ArrayListAppendBoolean(___94, ___2086)          ArrayListAppendTypedItem((___94),(___2086),___372)
 #define ArrayListAppendArbParam(___94, ___2086)         ArrayListAppendTypedItem((___94),(___2086),___90)
 #define ArrayListAppendUInt8Ptr(___94, ___2086)         ArrayListAppendTypedItem((___94),(___2086),uint8_t*)
 #define ArrayListAppendUInt16tPtr(___94, ___2086)       ArrayListAppendTypedItem((___94),(___2086),uint16_t*)
 #define ArrayListAppendUInt32tr(___94, ___2086)         ArrayListAppendTypedItem((___94),(___2086),uint32_t*)
 #define ArrayListAppendUInt64Ptr(___94, ___2086)        ArrayListAppendTypedItem((___94),(___2086),uint64_t*)
 #define ArrayListAppendInt64Ptr(___94, ___2086)         ArrayListAppendTypedItem((___94),(___2086),int64_t*)
 #define ArrayListAppendCharPtr(___94, ___2086)          ArrayListAppendTypedItem((___94),(___2086),char*)
 #define ArrayListAppendShortPtr(___94, ___2086)         ArrayListAppendTypedItem((___94),(___2086),short*)
 #define ArrayListAppendIntPtr(___94, ___2086)           ArrayListAppendTypedItem((___94),(___2086),int*)
 #define ArrayListAppendLongPtr(___94, ___2086)          ArrayListAppendTypedItem((___94),(___2086),long*)
 #define ArrayListAppendFloatPtr(___94, ___2086)         ArrayListAppendTypedItem((___94),(___2086),float*)
 #define ArrayListAppendDoublePtr(___94, ___2086)        ArrayListAppendTypedItem((___94),(___2086),double*)
 #define ArrayListAppendLgIndexPtr(___94, ___2086)       ArrayListAppendTypedItem((___94),(___2086),___2227*)
 #define ArrayListAppendEntIndexPtr(___94, ___2086)      ArrayListAppendTypedItem((___94),(___2086),___1172*)
 #define ArrayListAppendSmIntegerPtr(___94, ___2086)     ArrayListAppendTypedItem((___94),(___2086),int32_t*)
 #define ArrayListAppendBooleanPtr(___94, ___2086)       ArrayListAppendTypedItem((___94),(___2086),___372*)
 #define ArrayListAppendArbParamPtr(___94, ___2086)      ArrayListAppendTypedItem((___94),(___2086),___90*)
 #define ArrayListAppendVoidPtr(___94, ___2086)          ArrayListAppendTypedItem((___94),(___2086),void*)
 #define ArrayListAppendFunctionPtr(___94, ___2086)      ArrayListAppendTypedItem((___94),(___2086),(**)(void))
 #endif 
