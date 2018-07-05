function DataTable(data, titles, identifier) {

    this.data = data;
    this.titles = titles;
    this.identifier = identifier;

    this.initialize_data_table = function () {
        var dataTable = $(this.identifier).DataTable({
            "lengthMenu": [[50, 100, -1], [50, 100, "All"]],
            data: this.data,
            columns: this.titles,
            searching: true,
            autoWidth: false,
            dom: '<"table_buttons"B>l<"result_table"t><"table_pagination"p>',
            columnDefs: [
                {
                    targets: "_all",
                    className: 'dt-body-center'
                },
                {
                    targets: "_all",
                    className: 'dt-head-center'
                },
            ],
            buttons: [
                {
                    extend: 'csv',
                    text: 'Download table as csv'
                },
                {
                    extend: 'excel',
                    text: 'Download table as Excel (.xlsx)'
                }
            ],
        });
        return dataTable;
    };

    this.dataTable = this.initialize_data_table();

    this.addFiltering = function(bin_filter_identifier, suggestions, column) {
        $(bin_filter_identifier).puiautocomplete({
            completeSource: suggestions,
            multiple: true,
        });
        var that = this;
        $(bin_filter_identifier).on('focusin focusout keyup', function () {
            var all_items_htmlcoll = this.parentElement.parentElement.getElementsByTagName('li');
            var all_items = Array.prototype.slice.call(all_items_htmlcoll);
            //makes a collection of the raw text, removes empty entries and builds a regex string out of this collection
            var search_expr = all_items
                                .map(x => x.textContent)
                                .filter(x => x.length > 0)
                                .map(x => '^' + x + '$')
                                .join("|");
            that.dataTable
                .columns(column)
                .search(search_expr, true, false, true)
                .draw();
        });
    }

    this.add_colvis_filter = function(button_text, many_columns) {
        var that = this;
        if(many_columns) {
            var collectionLayout = "fixed four-column";
        }
        else {
            var collectionLayout = "fixed";
        }
        var colvisOptions = {
            extend: "colvis",
            collectionLayout: collectionLayout,
            text: button_text,
            prefixButtons: [
                {text : "Show all / Hide all", className: "show_hide_all_button",action: function ( e, dt, node, config ) {
                   var active = this.active();
                   this.active(!active);
                    this.columns().every(function() {
                        this.visible(active);
                    });
                }}
            ]};
        this.dataTable.button().add(0, colvisOptions);
    }
}


function performAjax(form_identifier, url, update_components) {
    console.log("serialized data: ");
    console.log($(form_identifier).serialize());
    $.ajaxSetup({
        headers: { "X-CSRFToken": Cookies.get('csrftoken')}
    });
    $.ajax({
        url: "update/",
        dataType: "json",
        data: $('#update_result_form').serialize(),
        success: function(result) {
            return result;
        },
        error: function() {
            console.error("Ajax request failed");
            console.trace();
            return null;
        }
    });
}

function DataTableData() {
    this.prediction_details_values;
    this.prediction_details_titles;
    this.bin_alias_list;
    this.model_list;
    this.prediction_values;
    this.prediction_titles;
    this.trait_counts_values;
    this.trait_counts_titles;
    this.bin_summary_values;
    this.bin_summary_titles;
}

function InitializeAllDataTables(dataTableData) {

    this.dataTableData = dataTableData;
    this.tables = {};
    /*this.table_prediction_details = null;
    this.table_prediction = null;
    this.table_trait_counts = null;
    this.table_bin_summary = null; */

    this.update = function(dataTableData) {
        this.dataTableData = dataTableData;
        for(var key in this.tables) {
            this.tables[key].dataTable.clear();
        }
        this.initialize();
        for(var key in this.tables) {
            this.tables[key].dataTable.draw()
        }
    };

    this.initialize = function()
    {
        try {
            this.tables['table_prediction_details'] = new DataTable(this.dataTableData.prediction_details_values, this.dataTableData.prediction_details_titles, "#table_prediction_details");
            this.tables['table_prediction_details'].addFiltering("#trait_prediction_details_bin_filter", this.dataTableData.bin_alias_list, 0);
            this.tables['table_prediction_details'].addFiltering('#trait_prediction_details_model_filter', this.dataTableData.model_list, 1)
            this.tables['table_prediction_details'].add_colvis_filter("Show/Hide columns", false);
        } catch (e) {
            console.error("table_prediction_details failed" + " " + e);
        }
        try {
            this.tables['table_prediction'] = new DataTable(this.dataTableData.prediction_values, this.dataTableData.prediction_titles, "#table_prediction");
            this.tables['table_prediction'].add_colvis_filter("Filter Model Columns", true);
        } catch (e) {
            console.error("table_prediction failed" + " " + e);
        }
        try {
            this.tables['table_trait_counts'] = new DataTable(this.dataTableData.trait_counts_values, this.dataTableData.trait_counts_titles, "#table_trait_counts");
            this.tables['table_trait_counts'].addFiltering("#trait_count_model_filter", this.dataTableData.model_list, 0)
        } catch (e) {
            console.error("table_trait_counts failed" + " " + e);
        }
        try {
            this.tables['table_bin_summary'] = new DataTable(this.dataTableData.bin_summary_values, this.dataTableData.bin_summary_titles, "#table_bin_summary");
            this.tables['table_bin_summary'].addFiltering("#trait_count_bin_filter", this.dataTableData.bin_alias_list, 0);
        } catch (e) {
            console.error("table_bin_summary failed" + " " + e);
        }
    }
}