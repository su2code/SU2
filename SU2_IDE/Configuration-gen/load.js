function load(){
	var resultJson = "";
	jQuery.ajaxSetup({async:false});

	$.get('option_structure.hpp', function(data) {

            var lines = data.split("\n");

			var start = false;
			var end = false;
			var blockStart = false;
			var blockEnd = false;
			var content = '\"FormBuilder\" : { \"fields\" : {';
            $.each(lines, function(n, elem) {
                if(start && !end){
                	if(blockStartFn(elem)) {
                		blockStart = true;
                		blockEnd = false;
                	}
                	if(blockStart && !blockEnd){
                		var groupName;
                		var listValue;
                		if(blockStartFn(elem)){
                			var res = elem.split(",");
                			groupName = res[2].trim();
                			groupName = groupName.substring(0, groupName.length - 1);

	                		content += '\"list_' + groupName + '\" : { \"label\" : \"' + groupName + '\", \"type\" : \"select\", \"options\" : {';
                		} else {
                			var res1 = elem.split(",");
                			listValue = res1[1].trim();
                			if(blockEndFn(elem)){
                				listValue = listValue.substring(0, listValue.length - 2);
                			} else {
                				listValue = listValue.substring(0, listValue.length - 1);
                		    }
                			content += '\"' + listValue + '\" : { \"label\" : \"' + listValue + '\"},';
                			if(blockEndFn(elem)){
                				content = content.substring(0, content.length-1);
	                			content += '}},';
                			}
                		}
                	}

                	if(blockStart == true){
                		if(blockEndFn(elem)){
                			blockStart = false;
                			blockEnd = true;
                		}
                	}

                }
                if(elem.indexOf("BEGIN_CONFIG_ENUMS") > -1){
                	start = true;
                }
                if(elem.indexOf("END_CONFIG_ENUMS") > -1){
                	end = true;
                }
            });
            content = content.substring(0, content.length-1);
			resultJson = '{' + content + '}}}';
        }, "text");
        return resultJson;
}

function blockStartFn(elem) {
	return elem.startsWith('static const map');
}

function blockEndFn(elem) {
	return elem.indexOf(';') > -1;
}
