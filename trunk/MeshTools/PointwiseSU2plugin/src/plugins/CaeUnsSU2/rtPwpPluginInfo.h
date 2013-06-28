/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

#ifndef _RTPWPPLUGININFO_H_
#define _RTPWPPLUGININFO_H_

/*! \cond */

    /* initialize the PWP_PLUGININFO data returned by
       PwpGetPluginInfo(PWP_PLUGININFO *pInfo)
    */
    VERSION_PWP_INIT,         // conforms to this PWP-API version
    VERSION_LIB_INIT,         // software library release version
    "Stanford University Aerospace Design Lab",        // company/author description
    "http://adl.stanford.edu/docs/display/SUSQUARED/Contact (susquared-dev@lists.stanford.edu)",  // support description (phone, web-link).
    "Copyright(c) 2012", // copyright description
    0,                        // number of APIs (auto-set at runtime)
    0,                        // default msg callback (auto-set at runtime)
    0,                        // spy msg callback (auto-set at runtime)

/*! \endcond */

/************************************************************************/
/*! \file
\brief Static Initialization Data for the PWP_PLUGININFO structure

initialize the PWP_PLUGININFO data 

The file \sf{%rtPwpPluginInfo.h} defines the static, compile-time
initialization of the PWP_PLUGININFO data struct returned by the PWP-API
function PwpGetPluginInfo(). If you want to see the implementation details,
look in the \sf{/shared/PWP/apiPWP.c} file.

The SDK file \sf{/shared/PWP/apiPWP.c} includes \sf{%rtPwpPluginInfo.h} as
shown below.
\par
\dontinclude apiPWP.c
\skip PWP_PLUGININFO info =
\until };

The format of \sf{%rtPwpPluginInfo.h} must be valid for the static
initialization of a C-struct. If you are not familiar with static
initialization, see the \ref example_cstruct_init page.

When copied from the \sf{src/plugins/templates/PWP/} folder to your plugins
project folder, \sf{%rtPwpPluginInfo.h} will contain example initilization
data. This example data must be edited to define the values appropriate for
your plugin's implementation.

\par Example PWP_PLUGININFO data

The example data from the SDK \sf{%rtPwpPluginInfo.h} template file is shown
below.
\par
\dontinclude rtPwpPluginInfo.h
\skip VERSION_PWP_INIT
\until spy msg callback

Please notice that the last three data members are initilized to 0. These
values are set automatically at runtime by PwpGetPluginInfo() as shown in the
code below.
\par
\dontinclude apiPWP.c
\skip info.apiCount
\until info.spyCB
*/

#endif /* _RTPWPPLUGININFO_H_ */
