angular.module('app', ['dynform'])
    .controller('AppCtrl', ['$scope', function ($scope) {
        $scope.urlFormData = {};   // JavaScript needs an object to put our form's models into.

        $scope.processForm = function () {
            alert ($scope.urlFormData.email);
        };
    }])
    .filter('pretty', function() {
        return function (input) {
            var temp;
            try {
                //temp = angular.fromJson(input);
                temp = input.email
            }
            catch (e) {
                temp = input;
            }

            return angular.toJson(temp, true);
        };
    });
