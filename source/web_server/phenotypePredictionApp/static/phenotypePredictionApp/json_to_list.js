
function json_to_datatable_input(json_obj) {
    var all_data = [];
    var titles = []

    var all_keys = Object.keys(json_obj[0].fields);
    for(key in all_keys) {
        let dicTmp = {title : key};
        titles.push(dicTmp);
    }

    for (i=0; i<json_obj.length; i++) {
        var sub_arr = Object.values(json_obj[i].fields);
        result_arr.push(sub_arr);
    }

    return [titles, all_data];
}