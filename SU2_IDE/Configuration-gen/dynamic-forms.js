/**
* DynamicForms - Build Forms in AngularJS From Nothing But JSON
* @version v0.0.4 - 2014-11-16
* @link http://github.com/danhunsaker/angular-dynamic-forms
* @license MIT, http://opensource.org/licenses/MIT
*/

/**
* Dynamically build an HTML form using a JSON array/object as a template.
*
* @todo Properly describe this directive.
* @param {mixed} [template] - The form template itself, as an array or object.
* @param {string} [templateUrl] - The URL to retrieve the form template from; template overrides.
* @param {Object} ngModel - An object in the current scope where the form data should be stored.
* @example <dynamic-form template-url="form-template.js" ng-model="formData"></dynamic-form>
*/
angular.module('dynform', [])
  .directive('dynamicForm', ['$q', '$parse', '$http', '$templateCache', '$compile', '$document', '$timeout', function ($q, $parse, $http, $templateCache, $compile, $document, $timeout) {
    var supported = {
        //  Text-based elements
        'text': {element: 'input', type: 'text', editable: true, textBased: true},
        'date': {element: 'input', type: 'date', editable: true, textBased: true},
        'datetime': {element: 'input', type: 'datetime', editable: true, textBased: true},
        'datetime-local': {element: 'input', type: 'datetime-local', editable: true, textBased: true},
        'email': {element: 'input', type: 'email', editable: true, textBased: true},
        'month': {element: 'input', type: 'month', editable: true, textBased: true},
        'number': {element: 'input', type: 'number', editable: true, textBased: true},
        'password': {element: 'input', type: 'password', editable: true, textBased: true},
        'search': {element: 'input', type: 'search', editable: true, textBased: true},
        'tel': {element: 'input', type: 'tel', editable: true, textBased: true},
        'textarea': {element: 'textarea', editable: true, textBased: true},
        'time': {element: 'input', type: 'time', editable: true, textBased: true},
        'url': {element: 'input', type: 'url', editable: true, textBased: true},
        'week': {element: 'input', type: 'week', editable: true, textBased: true},
        //  Specialized editables
        'checkbox': {element: 'input', type: 'checkbox', editable: true, textBased: false},
        'color': {element: 'input', type: 'color', editable: true, textBased: false},
        'file': {element: 'input', type: 'file', editable: true, textBased: false},
        'range': {element: 'input', type: 'range', editable: true, textBased: false},
        'select': {element: 'select', editable: true, textBased: false},
        //  Pseudo-non-editables (containers)
        'checklist': {element: 'div', editable: false, textBased: false},
        'fieldset': {element: 'fieldset', editable: false, textBased: false},
        'radio': {element: 'div', editable: false, textBased: false},
        //  Non-editables (mostly buttons)
        'button': {element: 'button', type: 'button', editable: false, textBased: false},
        'hidden': {element: 'input', type: 'hidden', editable: false, textBased: false},
        'image': {element: 'input', type: 'image', editable: false, textBased: false},
        'legend': {element: 'legend', editable: false, textBased: false},
        'reset': {element: 'button', type: 'reset', editable: false, textBased: false},
        'submit': {element: 'button', type: 'submit', editable: false, textBased: false}
      };
    
    return {
      restrict: 'E', // supports using directive as element only
      link: function ($scope, element, attrs) {
        //  Basic initialization
        var newElement = null,
          newChild = null,
          optGroups = {},
          cbAtt = '',
          foundOne = false,
          iterElem = element,
          model = null;
        
        //  Check that the required attributes are in place
        if (angular.isDefined(attrs.ngModel) && (angular.isDefined(attrs.template) || angular.isDefined(attrs.templateUrl)) && !element.hasClass('dynamic-form')) {
          model = $parse(attrs.ngModel)($scope);
          //  Grab the template. either from the template attribute, or from the URL in templateUrl
          (attrs.template ? $q.when($parse(attrs.template)($scope)) :
            $http.get(attrs.templateUrl, {cache: $templateCache}).then(function (result) {
              return result.data;
            })
          ).then(function (template) {
            var setProperty = function (obj, props, value, lastProp, buildParent) {
              props = props.split('.');
              lastProp = lastProp || props.pop();
              
              for (var i = 0; i < props.length; i++) {
                obj = obj[props[i]] = obj[props[i]] || {};
              }
              
              if (!buildParent) {
                obj[lastProp] = value;
              }
            },
            bracket = function (model, base) {
              props = model.split('.');
              return (base || props.shift()) + (props.length ? "['" + props.join("']['") + "']" : '');
            },
            buildFields = function (field, id) {
              if (String(id).charAt(0) == '$') {
                // Don't process keys added by Angular...  See GitHub Issue #29
                return;
              }
              else if (!angular.isDefined(supported[field.type]) || supported[field.type] === false) {
                //  Unsupported.  Create SPAN with field.label as contents
                newElement = angular.element('<span></span>');
                if (angular.isDefined(field.label)) {angular.element(newElement).html(field.label);}
                angular.forEach(field, function (val, attr) {
                  if (["label", "type"].indexOf(attr) > -1) {return;}
                  newElement.attr(attr, val);
                });
                this.append(newElement);
                newElement = null;
              }
              else {
                //  Supported.  Create element (or container) according to type
                if (!angular.isDefined(field.model)) {
                  field.model = id;
                }
                
                newElement = angular.element($document[0].createElement(supported[field.type].element));
                if (angular.isDefined(supported[field.type].type)) {
                  newElement.attr('type', supported[field.type].type);
                }
                
                //  Editable fields (those that can feed models)
                if (angular.isDefined(supported[field.type].editable) && supported[field.type].editable) {
                  newElement.attr('name', bracket(field.model));
                  newElement.attr('ng-model', bracket(field.model, attrs.ngModel));
                  // Build parent in case of a nested model
                  setProperty(model, field.model, {}, null, true);
                    
                  if (angular.isDefined(field.readonly)) {newElement.attr('ng-readonly', field.readonly);}
                  if (angular.isDefined(field.required)) {newElement.attr('ng-required', field.required);}
                  if (angular.isDefined(field.val)) {
                    setProperty(model, field.model, angular.copy(field.val));
                    newElement.attr('value', field.val);
                  }
                }
                
                //  Fields based on input type=text
                if (angular.isDefined(supported[field.type].textBased) && supported[field.type].textBased) {
                  if (angular.isDefined(field.minLength)) {newElement.attr('ng-minlength', field.minLength);}
                  if (angular.isDefined(field.maxLength)) {newElement.attr('ng-maxlength', field.maxLength);}
                  if (angular.isDefined(field.validate)) {newElement.attr('ng-pattern', field.validate);}
                  if (angular.isDefined(field.placeholder)) {newElement.attr('placeholder', field.placeholder);}
                }
                
                //  Special cases
                if (field.type === 'number' || field.type === 'range') {
                  if (angular.isDefined(field.minValue)) {newElement.attr('min', field.minValue);}
                  if (angular.isDefined(field.maxValue)) {newElement.attr('max', field.maxValue);}
                  if (angular.isDefined(field.step)) {newElement.attr('step', field.step);}
                }
                else if (['text', 'textarea'].indexOf(field.type) > -1) {
                  if (angular.isDefined(field.splitBy)) {newElement.attr('ng-list', field.splitBy);}
                }
                else if (field.type === 'checkbox') {
                  if (angular.isDefined(field.isOn)) {newElement.attr('ng-true-value', field.isOn);}
                  if (angular.isDefined(field.isOff)) {newElement.attr('ng-false-value', field.isOff);}
                  if (angular.isDefined(field.slaveTo)) {newElement.attr('ng-checked', field.slaveTo);}
                }
                else if (field.type === 'checklist') {
                  if (angular.isDefined(field.val)) {
                    setProperty(model, field.model, angular.copy(field.val));
                  }
                  if (angular.isDefined(field.options)) {
                    if (! (angular.isDefined(model[field.model]) && angular.isObject(model[field.model]))) {
                      setProperty(model, field.model, {});
                    }
                    angular.forEach(field.options, function (option, childId) {
                      newChild = angular.element('<input type="checkbox" />');
                      newChild.attr('name', bracket(field.model + '.' + childId));
                      newChild.attr('ng-model', bracket(field.model + "." + childId, attrs.ngModel));
                      if (angular.isDefined(option['class'])) {newChild.attr('ng-class', option['class']);}
                      if (angular.isDefined(field.disabled)) {newChild.attr('ng-disabled', field.disabled);}
                      if (angular.isDefined(field.readonly)) {newChild.attr('ng-readonly', field.readonly);}
                      if (angular.isDefined(field.required)) {newChild.attr('ng-required', field.required);}
                      if (angular.isDefined(field.callback)) {newChild.attr('ng-change', field.callback);}
                      if (angular.isDefined(option.isOn)) {newChild.attr('ng-true-value', option.isOn);}
                      if (angular.isDefined(option.isOff)) {newChild.attr('ng-false-value', option.isOff);}
                      if (angular.isDefined(option.slaveTo)) {newChild.attr('ng-checked', option.slaveTo);}
                      if (angular.isDefined(option.val)) {
                        setProperty(model, field.model, angular.copy(option.val), childId);
                        newChild.attr('value', option.val);
                      }
                      
                      if (angular.isDefined(option.label)) {
                          newChild = newChild.wrap('<label></label>').parent();
                          newChild.append(document.createTextNode(' ' + option.label));
                      }
                      newElement.append(newChild);
                    });
                  }
                }
                else if (field.type === 'radio') {
                  if (angular.isDefined(field.val)) {
                    setProperty(model, field.model, angular.copy(field.val));
                  }
                  if (angular.isDefined(field.values)) {
                    angular.forEach(field.values, function (label, val) {
                      newChild = angular.element('<input type="radio" />');
                      newChild.attr('name', bracket(field.model));
                      newChild.attr('ng-model', bracket(field.model, attrs.ngModel));
                      if (angular.isDefined(field['class'])) {newChild.attr('ng-class', field['class']);}
                      if (angular.isDefined(field.disabled)) {newChild.attr('ng-disabled', field.disabled);}
                      if (angular.isDefined(field.callback)) {newChild.attr('ng-change', field.callback);}
                      if (angular.isDefined(field.readonly)) {newChild.attr('ng-readonly', field.readonly);}
                      if (angular.isDefined(field.required)) {newChild.attr('ng-required', field.required);}
                      newChild.attr('value', val);
                      if (angular.isDefined(field.val) && field.val === val) {newChild.attr('checked', 'checked');}
                      
                      if (label) {
                          newChild = newChild.wrap('<label></label>').parent();
                          newChild.append(document.createTextNode(' ' + label));
                      }
                      newElement.append(newChild);
                    });
                  }
                }
                else if (field.type === 'select') {
                  if (angular.isDefined(field.multiple) && field.multiple !== false) {newElement.attr('multiple', 'multiple');}
                  if (angular.isDefined(field.empty) && field.empty !== false) {newElement.append(angular.element($document[0].createElement('option')).attr('value', '').html(field.empty));}
                  
                  if (angular.isDefined(field.autoOptions)) {
                    newElement.attr('ng-options', field.autoOptions);
                  }
                  else if (angular.isDefined(field.options)) {
                    angular.forEach(field.options, function (option, childId) {
                      newChild = angular.element($document[0].createElement('option'));
                      newChild.attr('value', childId);
                      if (angular.isDefined(option.disabled)) {newChild.attr('ng-disabled', option.disabled);}
                      if (angular.isDefined(option.slaveTo)) {newChild.attr('ng-selected', option.slaveTo);}
                      if (angular.isDefined(option.label)) {newChild.html(option.label);}
                      if (angular.isDefined(option.group)) {
                        if (!angular.isDefined(optGroups[option.group])) {
                          optGroups[option.group] = angular.element($document[0].createElement('optgroup'));
                          optGroups[option.group].attr('label', option.group);
                        }
                        optGroups[option.group].append(newChild);
                      }
                      else {
                        newElement.append(newChild);
                      }
                    });
                    
                    if (!angular.equals(optGroups, {})) {
                      angular.forEach(optGroups, function (optGroup) {
                        newElement.append(optGroup);
                      });
                      optGroups = {};
                    }
                  }
                }
                else if (field.type === 'image') {
                  if (angular.isDefined(field.label)) {newElement.attr('alt', field.label);}
                  if (angular.isDefined(field.source)) {newElement.attr('src', field.source);}
                }
                else if (field.type === 'hidden') {
                  newElement.attr('name', bracket(field.model));
                  newElement.attr('ng-model', bracket(field.model, attrs.ngModel));
                  if (angular.isDefined(field.val)) {
                    setProperty(model, field.model, angular.copy(field.val));
                    newElement.attr('value', field.val);
                  }
                }
                else if (field.type === 'file') {
                  if (angular.isDefined(field.multiple)) {
                    newElement.attr('multiple', field.multiple);
                  }
                }
                else if (field.type === 'fieldset') {
                  if (angular.isDefined(field.fields)) {
                    var workingElement = newElement;
                    angular.forEach(field.fields, buildFields, newElement);
                    newElement = workingElement;
                  }
                }
                
                //  Common attributes; radio already applied these...
                if (field.type !== "radio") {
                  if (angular.isDefined(field['class'])) {newElement.attr('ng-class', field['class']);}
                  //  ...and checklist has already applied these.
                  if (field.type !== "checklist") {
                    if (angular.isDefined(field.disabled)) {newElement.attr('ng-disabled', field.disabled);}
                    if (angular.isDefined(field.callback)) {
                      //  Some input types need listeners on click...
                      if (["button", "fieldset", "image", "legend", "reset", "submit"].indexOf(field.type) > -1) {
                        cbAtt = 'ng-click';
                      }
                      //  ...the rest on change.
                      else {
                        cbAtt = 'ng-change';
                      }
                      newElement.attr(cbAtt, field.callback);
                    }
                  }
                }
                
                //  If there's a label, add it.
                if (angular.isDefined(field.label)) {
                  //  Some elements have already applied their labels.
                  if (["image", "hidden"].indexOf(field.type) > -1) {
                    angular.noop();
                  }
                  //  Fieldset elements put their labels in legend child elements.
                  else if (["fieldset"].indexOf(field.type) > -1) {
                    newElement.prepend(angular.element($document[0].createElement('legend')).html(field.label));
                  }
                  //  Button elements get their labels from their contents.
                  else if (["button", "legend", "reset", "submit"].indexOf(field.type) > -1) {
                    newElement.html(field.label);
                  }
                  //  Everything else should be wrapped in a label tag.
                  else {
                    newElement = newElement.wrap('<label></label>').parent();
                    newElement.prepend(document.createTextNode(field.label + ' '));
                  }
                }
                
                // Arbitrary attributes
                if (angular.isDefined(field.attributes)) {
                  angular.forEach(field.attributes, function (val, attr) {
                    newElement.attr(attr, val);
                  });
                }
                
                // Add the element to the page
                this.append(newElement);
                newElement = null;
              }
            };
            
            angular.forEach(template, buildFields, element);
            
            //  Determine what tag name to use (ng-form if nested; form if outermost)
            while (!angular.equals(iterElem.parent(), $document) && !angular.equals(iterElem[0], $document[0].documentElement)) {
              if (['form','ngForm','dynamicForm'].indexOf(attrs.$normalize(angular.lowercase(iterElem.parent()[0].nodeName))) > -1) {
                foundOne = true;
                break;
              }
              iterElem = iterElem.parent();
            }
            if (foundOne) {
              newElement = angular.element($document[0].createElement('ng-form'));
            }
            else {
              newElement = angular.element("<form></form>");
            }
            
            //  Psuedo-transclusion
            angular.forEach(attrs.$attr, function(attName, attIndex) {
              newElement.attr(attName, attrs[attIndex]);
            });
            newElement.attr('model', attrs.ngModel);
            newElement.removeAttr('ng-model');
            angular.forEach(element[0].classList, function(clsName) {
              newElement[0].classList.add(clsName);
            });
            newElement.addClass('dynamic-form');
            newElement.append(element.contents());
            
            //  onReset logic
            newElement.data('$_cleanModel', angular.copy(model));
            newElement.bind('reset', function () {
              $timeout(function () {
                $scope.$broadcast('reset', arguments);
              }, 0);
            });
            $scope.$on('reset', function () {
              $scope.$apply(function () {
                $scope[attrs.ngModel] = {};
              });
              $scope.$apply(function () {
                $scope[attrs.ngModel] = angular.copy(newElement.data('$_cleanModel'));
              });
            });
            
            //  Compile and update DOM
            $compile(newElement)($scope);
            element.replaceWith(newElement);
          });
        }
      }
    };
  }])
  //  Not a fan of how Angular's ngList is implemented, so here's a better one (IMO).  It will ONLY
  //  apply to <dynamic-form> child elements, and replaces the ngList that ships with Angular.
  .directive('ngList', [function () {
    return {
      require: '?ngModel',
      link: function (scope, element, attr, ctrl) {
        var match = /\/(.*)\//.exec(element.attr(attr.$attr.ngList)),
          separator = match && new RegExp(match[1]) || element.attr(attr.$attr.ngList) || ',';
        
        if (element[0].form !== null && !angular.element(element[0].form).hasClass('dynamic-form')) {
          return;
        }
        
        ctrl.$parsers.splice(0, 1);
        ctrl.$formatters.splice(0, 1);
        
        ctrl.$parsers.push(function(viewValue) {
          var list = [];
          
          if (angular.isString(viewValue)) {
            //  Don't have Angular's trim() exposed, so let's simulate it:
            if (String.prototype.trim) {
              angular.forEach(viewValue.split(separator), function(value) {
                if (value) list.push(value.trim());
              });
            }
            else {
              angular.forEach(viewValue.split(separator), function(value) {
                if (value) list.push(value.replace(/^\s*/, '').replace(/\s*$/, ''));
              });
            }
          }
          
          return list;
        });

        ctrl.$formatters.push(function(val) {
          var joinBy = angular.isString(separator) && separator || ', ';
          
          if (angular.isArray(val)) {
            return val.join(joinBy);
          }
          
          return undefined;
        });
      }
    };
  }])
  //  Following code was adapted from http://odetocode.com/blogs/scott/archive/2013/07/05/a-file-input-directive-for-angularjs.aspx
  .directive('input', ['$parse', function ($parse) {
    return {
      restrict: 'E',
      require: '?ngModel',
      link: function (scope, element, attrs, ctrl) {
        if (!ctrl) {
          // Doesn't have an ng-model attribute; nothing to do here.
          return;
        }
        
        if (attrs.type === 'file') {
          var modelGet = $parse(attrs.ngModel),
            modelSet = modelGet.assign,
            onChange = $parse(attrs.onChange),
            updateModel = function () {
              scope.$apply(function () {
                modelSet(scope, element[0].files);
                onChange(scope);
              });                    
            };
          
          ctrl.$render = function () {
            element[0].files = this.$viewValue;
          };
          element.bind('change', updateModel);
        }
        else if (attrs.type === 'range') {
          ctrl.$parsers.push(function (val) {
            if (val) {
              return parseFloat(val);
            }
          });
        }
      }
    };
  }])
  //  Following code was adapted from http://odetocode.com/blogs/scott/archive/2013/07/03/building-a-filereader-service-for-angularjs-the-service.aspx
  .factory('fileReader', ['$q', function ($q) {
    var onLoad = function(reader, deferred, scope) {
        return function () {
          scope.$apply(function () {
            deferred.resolve(reader.result);
          });
        };
      },
      onError = function (reader, deferred, scope) {
        return function () {
          scope.$apply(function () {
            deferred.reject(reader.error);
          });
        };
      },
      onProgress = function(reader, scope) {
        return function (event) {
          scope.$broadcast('fileProgress',
            {
              total: event.total,
              loaded: event.loaded,
              status: reader.readyState
            });
        };
      },
      getReader = function(deferred, scope) {
        var reader = new FileReader();
        reader.onload = onLoad(reader, deferred, scope);
        reader.onerror = onError(reader, deferred, scope);
        reader.onprogress = onProgress(reader, scope);
        return reader;
      };

    return {
      readAsArrayBuffer: function (file, scope) {
        var deferred = $q.defer(),
          reader = getReader(deferred, scope);         
        reader.readAsArrayBuffer(file);
        return deferred.promise;
      },
      readAsBinaryString: function (file, scope) {
        var deferred = $q.defer(),
          reader = getReader(deferred, scope);         
        reader.readAsBinaryString(file);
        return deferred.promise;
      },
      readAsDataURL: function (file, scope) {
        var deferred = $q.defer(),
          reader = getReader(deferred, scope);         
        reader.readAsDataURL(file);
        return deferred.promise;
      },
      readAsText: function (file, scope) {
        var deferred = $q.defer(),
          reader = getReader(deferred, scope);         
        reader.readAsText(file);
        return deferred.promise;
      }
    };
  }]);

/*  End of dynamic-forms.js */
