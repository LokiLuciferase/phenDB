
function resultList_to_datatable_input(json_obj) {
    var all_data = [];
    var titles = [];

    var all_keys = Object.keys(json_obj[0].fields);
    for(var i=0; i<all_keys.length; i++) {
        var dicTmp = {title : all_keys[i]};
        titles.push(dicTmp);
    }

    for (var i=0; i<json_obj.length; i++) {
        var sub_arr = Object.values(json_obj[i].fields);
        all_data.push(sub_arr);
    }

    return [titles, all_data];
}

//WARNING: This code could break if the Django model is changed
function all_models_to_list(json_obj) {
    var model_names = [];
    var model_descriptions = [];
    for(var i=0; i<json_obj.length; i++) {
        var row = Object.values(json_obj[i].fields);
        if(model_names.indexOf(row[0]) != -1) continue;
        model_names.push(row[0]);
        model_descriptions.push(row[2]);
    }
    return [model_names, model_descriptions];
}

//TODO: put in other file
function models_to_infotext(titles, descriptions) {
    var infotable = "<table>";
    for(var i=0; i<titles.length; i++) {
        infotable += "<tr>";
        for(var a=0; a<2; a++) {
            infotable += "<td>";
            infotable += titles[i];
            infotable += "</td>";
            infotable += "<td>";
            infotable += descriptions[i];
            infotable += "</td>";
            infotable += "</tr>";
            ++i;
        }
    }
    infotable += "</table>";
    return infotable;
}