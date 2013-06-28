/****************************************************************************
 *
 * CAEP Plugin for Stanford University Unstructured
 *
 * Author - Trent Lukaczyk, Aerospace Design Lab
 *          Contact: susquared-dev@lists.stanford.edu
 *          March 2012
 *
 ***************************************************************************/

#include "apiCAEP.h"
#include "apiCAEPUtils.h"
#include "apiPWP.h"
#include "runtimeWrite.h"

//------------------------------------------------
static int
caeuDimensionToInt(CAEP_ENUM_DIMENSION dim)
// convert pointwise dimension enum to integer
{
    switch (dim) {
        case PWP_DIMENSION_2D: return 2;
        case PWP_DIMENSION_3D: return 3;
    }
    return 0;
}

//------------------------------------------------
static void
writeDimension (CAEP_RTITEM *pRti)
// write dimension info
{
    // comment header
    fputs("% \n% Problem dimension \n% \n", pRti->fp);
    // print number of dimensions
    fprintf( pRti->fp, "NDIME= %i\n",
             caeuDimensionToInt((pRti)->pWriteInfo->dimension) );
}

//------------------------------------------------
static int
elemType2Int(PWGM_ENUM_ELEMTYPE type)
// convert pointwise element type enum to integer
{
    int result = 0;
    switch (type) {
        case PWGM_ELEMTYPE_BAR:     result = 3;  break;
        case PWGM_ELEMTYPE_TRI:     result = 5;  break;
        case PWGM_ELEMTYPE_QUAD:    result = 9;  break;
        case PWGM_ELEMTYPE_TET:     result = 10; break;
        case PWGM_ELEMTYPE_HEX:     result = 12; break;
        case PWGM_ELEMTYPE_WEDGE:   result = 13; break;
        case PWGM_ELEMTYPE_PYRAMID: result = 14; break;
        default: result = 0; break;
    }
    return result;
}

//------------------------------------------------
static int
countElements(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model) 
// count global total number of elements
{

	  // total element counter
	  PWP_UINT32 eTotCnt = 0;

	  // check for valid data
	  if ( pRti && model ) {

		    // working variables
		    PWGM_ELEMCOUNTS eCounts ;								// elemet count bin, dummy variable
		    PWP_UINT32      iBlk     = 0;						// block index
		    PWGM_HBLOCK     hBlk    = PwModEnumBlocks(model, iBlk);	// block handle

		    // count all elements across all blocks
		    while ( PWGM_HBLOCK_ISVALID(hBlk) ) {
			      // add element count for block
			      eTotCnt += PwBlkElementCount(hBlk, &eCounts); 
			      // increment block
			      hBlk = PwModEnumBlocks(model, ++iBlk);         
		    }

	  } // if valid data

	  return eTotCnt;

}

//------------------------------------------------
static void
writeElemData(CAEP_RTITEM *pRti, PWGM_ELEMDATA *pElemData, PWP_UINT32 *eTotCnt)
// write element data line with or without global element index number
{
	  // check for valid data
    if ( pRti && pElemData ) {
		
		    // element local index
        PWP_UINT32 ii;
        
        // write element type
        fprintf( pRti->fp, "%2i ",
                 elemType2Int(pElemData->type) );

        // write vertex index numbers for current element
        for (ii=0; ii < pElemData->vertCnt; ++ii) {
            fprintf( pRti->fp, " %4lu",
                    (unsigned long)pElemData->index[ii] );
        }

		    // check for valid pointer before writing element count
		    if (0 != eTotCnt) {  
			      // write global element id 
			      // eTotCnt resides in writeBlocks()
			      fprintf( pRti->fp, " %4lu",
				            (unsigned long)((*eTotCnt)++) );
		    }

		    // print line return
        fputs("\n", pRti->fp);

    } // if valid data
	
}

//------------------------------------------------
static int
writeBlock(CAEP_RTITEM *pRti, PWGM_HBLOCK hBlk, PWP_UINT32 *eTotCnt)
// write a block
{
    // fail flag
	  int result = 0;

	  // check for valid data
    if (pRti && PWGM_HBLOCK_ISVALID(hBlk)) {

        PWGM_ELEMDATA   eData ;     // element data bin
		    PWGM_ELEMCOUNTS eCounts ;		// elemet count bin, dummy variable
        PWP_UINT32      eCnt = PwBlkElementCount(hBlk, &eCounts); // local element count;

		    // being major step, initialize substeps
        if ( caeuProgressBeginStep(pRti,eCnt) ) {

			      // reset counter
			      eCnt = 0;

			      // print elements
			      while ( PwElemDataMod(PwBlkEnumElements(hBlk, eCnt++), &eData) ) {
				        // write this element
				        writeElemData(pRti, &eData, eTotCnt);
				        // check for abort and incerement step
				        if (!caeuProgressIncr(pRti)) {
						        break;
				        }
			      } // while print elements

		    } // if beginStep

        // end progress step and check for abort
        result = caeuProgressEndStep(pRti);

    } // if valid data

	  // returns 1 if successful, 0 if not, give to while loop control
    return result;

}

//------------------------------------------------
static void
writeBlocks(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model)
// write multiple blocks
{
	    // check for user abort
      if (!pRti->opAborted) {

	        PWP_UINT32 cnt     = PwModBlockCount(model);    // block index, get number of blocks
	        PWP_UINT32 eTotCnt = countElements(pRti,model); // global element counter, get number of elements

	        // write comment header
          fputs("% \n% Inner element connectivity \n% \n", pRti->fp);

          // write number of elements
          fprintf( pRti->fp, "NELEM= %lu\n",
                  (unsigned long)eTotCnt );
		
	        // reset counts
          cnt     = 0;
	        eTotCnt = 0;

	        // write blocks, treat all blocks as one set of elements
          while ( writeBlock(pRti, PwModEnumBlocks(model, cnt++) , &eTotCnt ) ) {
		          // sneaky while loop ...
              // writeBlock() will break loop on user abort
          } 

      } // if aborted

}

//------------------------------------------------
static int
writeVertex(CAEP_RTITEM *pRti, PWGM_HVERTEX vertex)
// write a vertex
{
	  // fail flag
    int result = 0;

	  // check for valid data
    if (pRti && PWGM_HVERTEX_ISVALID(vertex)) {

        // print setup
        int prec = CAEPU_RT_PREC_SINGLE(pRti) ? 8 : 16; // print precision
        int spac = prec + 8;							// print spacing

        // get vertex data
        PWGM_VERTDATA VertData;
        PwVertDataMod(vertex, &VertData);  // pointer assigment...

        // print vertex data, 3D
        if ( CAEPU_RT_DIM_3D(pRti) ) {  // check dimension...
            fprintf( pRti->fp,
                     "%#*.*g %#*.*g %#*.*g %4lu\n",
                     spac, prec, VertData.x, 
                     spac, prec, VertData.y,
                     spac, prec, VertData.z, 
                    (unsigned long)VertData.i );
        } // if 3D print

        // print vertex data, 2D
        else {
            fprintf( pRti->fp,
                     "%#*.*g %#*.*g %4lu\n",
                     spac, prec, VertData.x, 
                     spac, prec, VertData.y,
                    (unsigned long)VertData.i );
        } // if 2D print

	      // check for success
        result = !pRti->opAborted;

    } // if valid data

	  // returns 1 if successful, 0 if not, give to while loop control
    return result;

}

//------------------------------------------------
static void
writeVertices(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model)
// write all vertex information
{
	  // check for user abort
    if (!pRti->opAborted) {

        // count number of verticies
        PWP_UINT32 cnt = PwModVertexCount(model);

        // being major step, initialize substeps
        if (caeuProgressBeginStep(pRti, cnt)) {

            // comment header
            fputs("% \n% Node coordinates \n% \n", pRti->fp);

            // number of points
            fprintf( pRti->fp, "NPOIN= %lu\n",
                    (unsigned long)cnt );

            // reset counter
            cnt=0;

	          // write verticies
            while ( writeVertex(pRti, PwModEnumVertices(model, cnt++)) ) {
                // sneaky while loop ...
                // writeVertex() will break loop on user abort
            }
			
	          // end steps
	          caeuProgressEndStep(pRti);

        } // if begin step

    } // if !abort

}

//------------------------------------------------
static void
writeCondData(CAEP_RTITEM *pRti, PWGM_CONDDATA *pCondData)
// write Condition Data (SU2 Marker Tag Name)
{
    // check for valid data
    if (pRti && pCondData) {
        // print Marker NameTag
        fprintf( pRti->fp, "MARKER_TAG= %s\n", 
				         pCondData->name );
    }
}

//------------------------------------------------
static int
writeDomain(CAEP_RTITEM *pRti, PWGM_HDOMAIN hDom)
// print a Domain (SU2 Marker)
{
    // fail flag
	  int result = 0;

	  // check for valid data
    if ( pRti && PWGM_HDOMAIN_ISVALID(hDom) ) {

        PWGM_ELEMDATA eData;									// element data bin
        PWGM_ELEMCOUNTS elemCnts;								// element count data bin
        PWP_UINT32 eCnt = PwDomElementCount(hDom, &elemCnts);   // get number of elements for marker
        PWGM_CONDDATA CondData;									// condition data (names and stuff)

        // Check if marker exists for this domain
        if ( PwDomCondition(hDom, &CondData) ) { 

		        // write marker top level info
            writeCondData(pRti, &CondData);

		        // initialize progress step
		        if (caeuProgressBeginStep(pRti, eCnt)) { 
            
			          // Write Number of Marker Elements
			          fprintf( pRti->fp, "MARKER_ELEMS= %lu\n",
				                (unsigned long)eCnt );

			          // reset counter
			          eCnt = 0;

			          // increment elements
			          while ( PwElemDataMod(PwDomEnumElements(hDom, eCnt++), &eData) ) {

				            // write elements (no index)
				            writeElemData(pRti, &eData, 0);

				            // increment step
				            if (!caeuProgressIncr(pRti)) {
					              break;
				            }

			          } // while increment elements

		        } // if initialize step

            // end progress step and check for abort
            result = caeuProgressEndStep(pRti);

        } // if domain exists

    } // if valid data

    return result;

}

//------------------------------------------------
static void
writeDomains(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model)
// write multiple Domains (SU2 Markers)
{
	  // check for user abort
    if (!pRti->opAborted) {

		    // marker counter
        PWP_UINT32 cnt = PwModDomainCount(model);
        
        // comment header
        fputs("% \n% Boundary elements \n% \n", pRti->fp);

        // number of markers
        fprintf( pRti->fp, "NMARK= %lu\n", 
                (unsigned long)cnt );

		    // reset counter
        cnt = 0;

		    // write markers
        while (writeDomain(pRti, PwModEnumDomains(model, cnt++))) {
            // sneaky while loop...
            // writeDomain() will break loop on user abort
        }

    } // if !abort

}

/**************************************/
// this is the main function
// it is called when you hit Save CAE in pointwise
PWP_BOOL
runtimeWrite( CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model,
              const CAEP_WRITEINFO *pWriteInfo )
{
	  // success flag
    PWP_BOOL result = PWP_FALSE;

    // check for valid data
    if (pRti && model) {
      
        // step information
        PWP_UINT32 nBlk     = PwModBlockCount(model);       // number of blocks
        PWP_UINT32 nDom     = PwModDomainCount(model);      // number of domains
        PWP_UINT32 stpMajor = 1 + nDom + nBlk;              // number of Major steps

        // initialize major steps
        if (caeuProgressInit(pRti, stpMajor)) {

            // Write Dimension
            // 0 major steps
            writeDimension (pRti);

            // Write Elements
            // nBlk major steps
            writeBlocks(pRti, model);

            // Write Points
            // 1 major steps
			      writeVertices(pRti, model);

            // Write Domains (SU2 Markers)
            // nDom major steps
            writeDomains(pRti, model);

            // check for fail
			      result = !pRti->opAborted;

        } // if initialize major steps

		    // Finish Major Steps
		    caeuProgressEnd(pRti, result);

    } // if valid data

    return result;

}
