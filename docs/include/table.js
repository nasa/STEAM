$(function (){
/* var createList = function(selector){
   
       var ul = $('<ul>');
       var selected = $(selector);
   
       if (selected.length === 0){
           return;
       }
   
       selected.clone().each(function (i,e){
   
           var p = $(e).children('dt').children('.descclassname');
           var n = $(e).children('dt').children('.descname');
           var l = $(e).children('dt').children('.headerlink');
           var d = $(e).children('dd').children('p').text();

           var a = $('<a>');
           a.attr('href',l.attr('href')).attr('title', 'Link to this definition');
   
           a.append(p).append(n);
   
           var entry = $('<li>').append(a).append(d);
           ul.append(entry);
       });
       return ul;
   }
   
   
   var c = $('<div style="float:left; width: 100%;">');
   
   var ul0 = c.clone().append($('.submodule-index'))
   
   customIndex = $('.custom-index');
   customIndex.empty();
   customIndex.append(ul0);
   
   var x = [];
   x.push(['Classes','dl.class']);
   x.push(['Functions','dl.function']);
   x.push(['Variables','dl.data']);
   x.push(['Exceptions','dl.exception']);
   
   x.forEach(function (e){
       var l = createList(e[1]);
       if (l) {        
           var ul = c.clone()
               .append('<p class="rubric">'+e[0]+'</p>')
               .append(l);
       }
       customIndex.append(ul);
   });
*/


   var createTable = function(selector){

       var table    = $("<table>")
         .attr('border','1')
         .addClass("longtable docutils")
         .append(`
<colgroup>
<col width="50%" />
<col width="50%" />
</colgroup>
   `);
// <tbody valign="top">

       var selected = $(selector);
   
       if (selected.length === 0){
           return;
       }
   
       selected.clone().each(function (i,e){
          var row      = $('<tr>');

          var linkCell = $('<td>');
          var textCell = $('<td>');
   
           var p = $(e).children('dt').children('.descclassname');
           var n = $(e).children('dt').children('.descname');
           var l = $(e).children('dt').children('.headerlink');
           var d = $(e).children('dd').children('p').text();

           var a = $('<a>');
           a.attr('href',l.attr('href')).attr('title', 'Link to this definition');
   
           a.append(p).append(n);
           linkCell.append(a)

           textCell.text(d);
   
           row.append(linkCell).append(textCell)
           
           table.append(row)
       });
       return table;
   }

   var c = $('<div style="float:left; width: 100%;">');

   customIndex = $('.custom-index');
   customIndex.empty();
   
   var x = [];
   x.push(['Classes','dl.class']);
   x.push(['Functions','dl.function']);
   x.push(['Variables','dl.data']);
   x.push(['Exceptions','dl.exception']);
   
   x.forEach(function (e){

       var out = createTable(e[1]);
       if (out) {
           var section = c.clone()
               .append('<p class="rubric">'+e[0]+'</p>')
               .append(out);
       }
          customIndex.append(section);

   });
   
});
