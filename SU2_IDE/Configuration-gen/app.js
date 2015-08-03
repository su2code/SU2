angular.module('app', ['dynform'])
    .controller('AppCtrl', ['$scope', '$http', function ($scope, $http) {

    $scope.urlFormData = {};   // JavaScript needs an object to put our form's models into.
	$scope.stdFormTemplate = load1();

	$scope.status={};
	$scope.status.message="";

	function load1(){
		var returnVal = '';
		jQuery.ajaxSetup({async:false});
		$.get('option_structure.hpp', function(data) {
			returnVal = {
      			"fieldset": {
        			"type": "fieldset",
        			"label": "fieldset",
        			"fields": {
        				"textIn": {
            				"type": "text",
            				"label": "textIn",
            				"placeholder": "text"
          				}
        			}
     			}
    		};
		}, 'text');
		return returnVal;
    }

        $scope.processForm = function () {
            alert ($scope.urlFormData.MACH_NUMBER);
        };
    }])
    .filter('nl2br', function() {
        return function (input) {
            var temp;
            try {
                //temp = angular.fromJson(input);
                temp = '<br/>' +
                        'MATH_PROBLEM=' + input.MATH_PROBLEM+ '<br/>' +
                        'PHYSICAL_PROBLEM=' +input.PHYSICAL_PROBLEM+ '<br />' +
'RESTART_SOL=' +input.RESTART_SOL+ '<br />' +
'AoA=' +input.AoA+ '<br />' +
'FREESTREAM_PRESSURE=' +input.FREESTREAM_PRESSURE+ '<br />' +
'FREESTREAM_TEMPERATURE=' +input.FREESTREAM_TEMP+ '<br />' +
'REF_ORIGIN_MOMENT_X=' +input.REF_ORIGIN_MOMENT_X+ '<br />' +
'REF_ORIGIN_MOMENT_Y=' +input.REF_ORIGIN_MOMENT_Y+ '<br />' +
'REF_ORIGIN_MOMENT_Z=' +input.REF_ORIGIN_MOMENT_Z+ '<br />' +
'REF_LENGTH_MOMENT=' +input.REF_LENGTH_MOMENT+ '<br />' +
'REF_AREA=' +input.REF_AREA+ '<br />' +
'REF_ELEM_LENGTH=' +input.REF_ELEM_LENGTH+ '<br />' +
'MARKER_EULER=' +input.MARKER_EULER+ '<br />' +
'MARKER_INLET=' +input.MARKER_INLET+ '<br />' +
'MARKER_OUTLET=' +input.MARKER_OUTLET+ '<br />' +
'MARKER_PLOTTING=' +input.MARKER_PLOTTING+ '<br />' +
'MARKER_MONITORING=' +input.MARKER_MONITORING+ '<br />' +
'NUM_METHOD_GRAD=' +input.NUM_METHOD_GRAD+ '<br />' +
'CFL_NUMBER=' +input.CFL_NUMBER+ '<br />' +
'CFL_ADAPT=' +input.CFL_ADAPT+ '<br />' +
'CFL_ADAPT_PARAM=' +input.CFL_ADAPT_PARAM+ '<br />' +
'RK_ALPHA_COEFF=' +input.RK_ALPHA_COEFF+ '<br />' +
'EXT_ITER=' +input.EXT_ITER+ '<br />' +
'LINEAR_SOLVER=' +input.LINEAR_SOLVER+ '<br />' +
'LINEAR_SOLVER_PREC=' +input.LINEAR_SOLVER_PREC+ '<br />' +
'LINEAR_SOLVER_ERROR=' +input.LINEAR_SOLVER_ERROR+ '<br />' +
'LINEAR_SOLVER_ITER=' +input.LINEAR_SOLVER_ITER+ '<br />' +
'MGLEVEL=' +input.MGLEVEL+ '<br />' +
'MGCYCLE=' +input.MGCYCLE+ '<br />' +
'MG_PRE_SMOOTH=' +input.MG_PRE_SMOOTH+ '<br />' +
'MG_POST_SMOOTH=' +input.MG_POST_SMOOTH+ '<br />' +
'MG_CORRECTION_SMOOTH=' +input.MG_CORRECTION_SMOOTH+ '<br />' +
'MG_DAMP_PROLONGATION=' +input.MG_DAMP_PROLONGATION+ '<br />' +
'CONV_NUM_METHOD_FLOW=' +input.CONV_NUM_METHOD_FLOW+ '<br />' +
'SPATIAL_ORDER_FLOW=' +input.SPATIAL_ORDER_FLOW+ '<br />' +
'SLOPE_LIMITER_FLOW=' +input.SLOPE_LIMITER_FLOW+ '<br />' +
'LIMITER_COEFF=' +input.LIMITER_COEFF+ '<br />' +
'AD_COEFF_FLOW=' +input.AD_COEFF_FLOW+ '<br />' +
'TIME_DISCRE_FLOW=' +input.TIME_DISCRE_FLOW+ '<br />' +
'CONV_CRITERIA=' +input.CONV_CRITERIA+ '<br />' +
'RESIDUAL_REDUCTION=' +input.RESIDUAL_REDUCTION+ '<br />' +
'RESIDUAL_MINVAL=' +input.RESIDUAL_MINVAL+ '<br />' +
'STARTCONV_ITER=' +input.STARTCONV_ITER+ '<br />' +
'CAUCHY_ELEMS=' +input.CAUCHY_ELEMS+ '<br />' +
'CAUCHY_EPS=' +input.CAUCHY_EPS+ '<br />' +
'CAUCHY_FUNC_FLOW=' +input.CAUCHY_FUNC_FLOW+ '<br />' +
'MESH_FILENAME=' +input.MESH_FILENAME+ '<br />' +
'MESH_FORMAT=' +input.MESH_FORMAT+ '<br />' +
'MESH_OUT_FILENAME=' +input.MESH_OUT_FILENAME+ '<br />' +
'SOLUTION_FLOW_FILENAME=' +input.SOLUTION_FLOW_FILENAME+ '<br />' +
'SOLUTION_LIN_FILENAME=' +input.SOLUTION_LIN_FILENAME+ '<br />' +
'SOLUTION_ADJ_FILENAME=' +input.SOLUTION_ADJ_FILENAME+ '<br />' +
'OUTPUT_FORMAT=' +input.OUTPUT_FORMAT+ '<br />' +
'CONV_FILENAME=' +input.CONV_FILENAME+ '<br />' +
'RESTART_FLOW_FILENAME=' +input.RESTART_FLOW_FILENAME+ '<br />' +
'RESTART_ADJ_FILENAME=' +input.RESTART_ADJ_FILENAME+ '<br />' +
'RESTART_LIN_FILENAME=' +input.RESTART_LIN_FILENAME+ '<br />' +
'VOLUME_FLOW_FILENAME=' +input.VOLUME_FLOW_FILENAME+ '<br />' +
'VOLUME_ADJ_FILENAME=' +input.VOLUME_ADJ_FILENAME+ '<br />' +
'VOLUME_LIN_FILENAME=' +input.VOLUME_LIN_FILENAME+ '<br />' +
'GRAD_OBJFUNC_FILENAME=' +input.GRAD_OBJFUNC_FILENAME+ '<br />' +
'SURFACE_FLOW_FILENAME=' +input.SURFACE_FLOW_FILENAME+ '<br />' +
'SURFACE_ADJ_FILENAME=' +input.SURFACE_ADJ_FILENAME+ '<br />' +
'SURFACE_LIN_FILENAME=' +input.SURFACE_LIN_FILENAME+ '<br />' +
'WRT_SOL_FREQ=' +input.WRT_SOL_FREQ+ '<br />' +
'WRT_CON_FREQ=' +input.WRT_CON_FREQ+
                		'<br/>';
            }
            catch (e) {
                temp = input;
            }

            return angular.toJson(temp, true);
        };
    });
