angular.module('app', ['dynform'])
    .controller('AppCtrl', ['$scope', '$http', function ($scope, $http) {

    $scope.urlFormData = {};   // JavaScript needs an object to put our form's models into.
	$scope.stdFormTemplate = load();

	$scope.status={};
	$scope.status.message="";

        $scope.processForm = function () {
            alert ($scope.urlFormData.MACH_NUMBER);
        };
    }])
    .filter('nl2br', function() {
        return function (input) {
            var temp;
            try {
                temp = '<br/>'

				$.each(input, function(k, v) {
					//display the key and value pair
					temp += k + '=' + v + '<br/>';
                });
                temp += '<br/>';

            }
            catch (e) {
                temp = input;
            }

            return angular.toJson(temp, true);
        };
    });
