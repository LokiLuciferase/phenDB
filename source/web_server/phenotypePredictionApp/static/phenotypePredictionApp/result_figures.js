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


function performAjax(form_id, url, update_components) {
    $.ajax({
        type: "POST",
        url: "update/",
        method: "POST",
        data: {"test_val" : "5"},
        contentType: "application/json",
        success: function(result) {
            console.log("ajax result: " + result);
        },
    });
}