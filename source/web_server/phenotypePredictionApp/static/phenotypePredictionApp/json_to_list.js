
function json_to_list(json_obj) {
    var result_arr = [];
    for (single_entry in json_obj) {
        var sub_arr = [];
        for (fieldKey in single_entry.fields) {
            console.log("single entry: " + single_entry.fields[fieldKey]);
            sub_arr.push(single_entry.fields[fieldKey]);
        }
        result_arr.push(sub_arr);
    }
    return result_arr;
}