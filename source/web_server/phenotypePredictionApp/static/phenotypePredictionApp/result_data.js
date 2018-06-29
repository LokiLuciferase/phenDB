var ResultData = {

    __data: null,

    initialize : function(data) {
        this.__data = data;
    },

    resultList_to_datatable_input : function() {
        var displayed_fields = ["bin", "model", "prediction", "confidence", "accuracy"];

        var returnval_data = [];
        var returnval_titles = [];


        var all_keys = Object.keys(json_obj[0].fields);
        titles = convertTitles(all_keys);

        for (var i=0; i<json_obj.length; i++) {
            var sub_arr = Object.values(json_obj[i].fields);
            all_data.push(sub_arr);
        }
        __roundNumbers(all_data, [3,4], 2);
        return [titles, all_data];
    },



}