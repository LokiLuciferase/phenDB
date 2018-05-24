
function json_to_list(json_obj) {
    var result_arr = [];
    var result_keys = Object.keys(json_obj[0].fields);
    for (i=0; i<json_obj.length; i++) {
        var sub_arr = Object.values(json_obj[i].fields);
        result_arr.push(sub_arr);
    }
    return [result_keys, result_arr];
}