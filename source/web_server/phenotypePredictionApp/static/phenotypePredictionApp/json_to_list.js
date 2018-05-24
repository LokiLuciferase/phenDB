
function json_to_list(json_obj) {
    var result_arr = [];
    for (i=0; i<json_obj.length; i++) {
        var sub_arr = Object.values(json_obj[i].fields);
        result_arr.push(sub_arr);
    }
    return result_arr;
}