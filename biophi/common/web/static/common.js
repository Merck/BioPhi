// Listen for "Load example" clicks


document.addEventListener('click', function (event) {
	if (!event.target.matches('a[data-inject-example]')) return;
	event.preventDefault();
    var inputSelector = event.target.getAttribute('href');
    var exampleText = event.target.getAttribute('data-inject-example');
    var inputElement = document.querySelector(inputSelector);
    inputElement.value = exampleText
}, false);

[].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]')).map(function (elem) {
  var args = {animation: false}
  if (elem.hasAttribute('data-tooltip-classes')) {
    args['customClass'] = elem.getAttribute('data-tooltip-classes');
  }
  return new bootstrap.Tooltip(elem, args)
})

document.querySelectorAll('tr[data-href]').forEach(elem => elem.addEventListener('click', function (event) {
    window.location = elem.getAttribute('data-href');
}))

document.querySelectorAll('[data-radio-set]').forEach(elem => elem.addEventListener('change', function (event) {
    var target = document.querySelector(elem.getAttribute('data-radio-set'));
    target.value = elem.getAttribute('data-value')
}))

document.querySelectorAll('[data-radio-enable]').forEach(elem => elem.addEventListener('change', function (event) {
    var target = document.querySelector(elem.getAttribute('data-radio-enable'));
    target.disabled = false
}))

document.querySelectorAll('[data-radio-disable]').forEach(elem => elem.addEventListener('change', function (event) {
    var target = document.querySelector(elem.getAttribute('data-radio-disable'));
    target.disabled = true
}))

document.querySelectorAll('[data-radio-tab]').forEach(elem => elem.addEventListener('change', function (event) {
    var target = document.querySelector(elem.getAttribute('data-radio-tab'));
    Array.from(target.parentNode.children).forEach(elem => elem.classList.remove('active'))
    target.classList.add('active')
}))

window.addEventListener("load", function(){
  document.querySelectorAll('[data-radio-set],[data-radio-enable],[data-radio-disable],[data-radio-tab]').forEach(function(elem) {
    if (elem.checked) {
      elem.dispatchEvent(new Event('change'))
    }
  })
}, false);

function isNumeric(str) {
  if (typeof str != "string") return false
  str = str.replace('%','')
  return !isNaN(str) && !isNaN(parseFloat(str))
}

document.querySelectorAll('th.sortable').forEach(elem => elem.addEventListener('click', sortTable))

// Taken from https://www.w3schools.com/howto/howto_js_sort_table.asp
function sortTable() {
  var th = this
  var tr = th.parentNode
  var rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
  var colNumber = Array.prototype.indexOf.call(tr.children, th);
  var table = th.closest('table');
  switching = true;
  // Set the sorting direction to ascending:
  dir = "asc";
  /* Make a loop that will continue until
  no switching has been done: */
  while (switching) {
    // Start by saying: no switching is done:
    switching = false;
    rows = table.rows;
    /* Loop through all table rows (except the
    first, which contains table headers): */
    for (i = 1; i < (rows.length - 1); i++) {
      // Start by saying there should be no switching:
      shouldSwitch = false;
      /* Get the two elements you want to compare,
      one from current row and one from the next: */
      var cells = rows[i].getElementsByTagName("td")
      /* Skip extra header rows */
      if (!cells.length) {
        continue
      }
      x = cells[colNumber];
      y = rows[i + 1].getElementsByTagName("td")[colNumber];
      /* Check if the two rows should switch place,
      based on the direction, asc or desc: */
      x_val = x.hasAttribute('data-value') ? x.getAttribute('data-value') : x.innerHTML.toLowerCase()
      y_val = y.hasAttribute('data-value') ? y.getAttribute('data-value') : y.innerHTML.toLowerCase()
      if (isNumeric(x_val) && isNumeric(y_val)) {
        x_val = parseFloat(x_val)
        y_val = parseFloat(y_val)
      }
      if (dir == "asc") {
        if (x_val > y_val) {
          // If so, mark as a switch and break the loop:
          shouldSwitch = true;
          break;
        }
      } else if (dir == "desc") {
        if (x_val < y_val) {
          // If so, mark as a switch and break the loop:
          shouldSwitch = true;
          break;
        }
      }
    }
    if (shouldSwitch) {
      /* If a switch has been marked, make the switch
      and mark that a switch has been done: */
      rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
      switching = true;
      // Each time a switch is done, increase this count by 1:
      switchcount ++;
    } else {
      /* If no switching has been done AND the direction is "asc",
      set the direction to "desc" and run the while loop again. */
      if (switchcount == 0 && dir == "asc") {
        dir = "desc";
        switching = true;
      }
    }
  }
  tr.querySelectorAll('th').forEach(elem => elem.classList.remove('sorted-asc', 'sorted-desc'))
  th.classList.add('sorted-'+dir)
}