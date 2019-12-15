function update_table(sid, what, x)
{
    request = new XMLHttpRequest();
    addr ='/update/' + sid + '/' + what + '/' + x + '/'
    var inputs = document.getElementsByTagName('input');
    sep = '?'
    for (var i = 0; i < inputs.length; i++) {
        addr += sep + inputs[i].name + '=' + inputs[i].value;
        sep = ',';
    }
    request.open('GET', addr, true);
    request.onload = function() {
        data = request.responseText;
        table = document.getElementById('database1')
        table.innerHTML = data;
    }
    request.send();
}
