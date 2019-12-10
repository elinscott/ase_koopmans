function update_table(sid, what, x)
{
    request = new XMLHttpRequest();
    q = document.getElementById('formula-result').value
    request.open('GET',
                 '/update/' + sid + '/' + what + '/' + x + '/?query=' + q,
                 true);
    request.onload = function() {
        data = request.responseText;
        table = document.getElementById('database1')
        table.innerHTML = data;
    }
    request.send();
}
