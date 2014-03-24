/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
#ifndef _W__BASE_h
#define _W__BASE_h

#if defined TECPLOTKERNEL

#include "C__BASE.h"

typedef ULONG_PTR Widget;
#define NULLWIDGET ((unsigned long)0)

extern BOOL MSWinIsLargeFont;
extern CWnd *(TPClassArray[TPNUMCLASSES]);

typedef struct
{
    CRect Rect;
    BOOL  RectSet;
    CRect MinRect;    // Only used for resizable dialogs
    BOOL  MinRectSet; // Only used for resizable dialogs
} TPDialogRect_s;

extern TPDialogRect_s TPDialogRect[TPNUMCLASSES];

#define WindowIsResizable(wnd) ((wnd)->GetStyle() & WS_THICKFRAME)

/*
 * NOTE: C = class number in TPClassArray (see C__BASE.h)
 *           0=no class
 *           1=TPMAINWND
 *       O = offset within container
 *           0 means the container itself
 *           any other number is the offset within control (starting at 1)
 *       I = control number/id
 *           0 = class itself
 *           for bitmap button groups, this is the container id
 *           for all other controls, this is the control id
 * WidgetValue 0 is the NULL widget.
 */

//#define MAKEWIDGET(C,I,O) ((Widget)((((C)&0xFFFF)<<16)|(((O)&0xFF)<<8)|((I)&0xFF)))

#define CLASS_NUM_BITS   10
#define ITEM_ID_BITS     16 /* must be at least 16-bits to fit in the workarea */
#define ITEM_OFFSET_BITS  6

/* The FILTER is also the maximum value. */
#define CLASS_NUM_FILTER   ((1<<CLASS_NUM_BITS)-1)
#define ITEM_ID_FILTER     ((1<<ITEM_ID_BITS)-1)
#define ITEM_OFFSET_FILTER ((1<<ITEM_OFFSET_BITS)-1)

#define TP_CLASS_NO_CLASS        0
#define TP_ITEM_CLASS_ITSELF     0
#define TP_OFFSET_ITEM_ITSELF    0

#define TP_OFFSET_MENU_CASCADE   (ITEM_OFFSET_FILTER-1)
#define TP_OFFSET_MENU_OPTION    ITEM_OFFSET_FILTER

/* MAKEWIDGET must be a define and not an inline function so widget's can be used in case statements. */
#define MAKEWIDGET(ClassNum, ItemID, ItemOffset) ((Widget)( ( ( ((ClassNum)<<ITEM_ID_BITS) | (ItemID) ) << ITEM_OFFSET_BITS ) | (ItemOffset) ))

inline INT_PTR GetClassNum(Widget W)
{
    return (W >> ITEM_ID_BITS) >> ITEM_OFFSET_BITS;
}
inline INT_PTR GetItemID(Widget W)
{
    return (W >> ITEM_OFFSET_BITS) & ITEM_ID_FILTER;
}
inline INT_PTR GetItemOffset(Widget W)
{
    return W & ITEM_OFFSET_FILTER;
}

/*
 * We cannot use Boolean_t here because it is not defined yet.
 */
extern CWnd  *GetWindowFromWidget(Widget);
extern CWnd  *GetParentWndFromWidget(Widget);
extern Widget GetWidgetFromWindow(CWnd *);
extern BOOL   WidgetExists(Widget W);
extern BOOL   TPEnableDialog(CWnd *wnd,
                             BOOL  bEnable);
extern void   TPEnableAllDialogs(BOOL bEnable);
extern TPDialogRect_s* GetDialogRectFromWidget(Widget);
extern TPDialogRect_s* GetDialogRectFromWindow(CWnd* pWnd);


#endif /* TECPLOTKERNEL */

/*
 * This is the start of the translated string ID's
 * in the resource file
 */
#define TP_TRANSLATED_STRING_START_ID 20000

#endif
