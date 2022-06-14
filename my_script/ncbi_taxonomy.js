cars = document.querySelector("#dummybodyid > table:nth-child(7) > tbody").childNodes
var i = 0;
temp = cars[i].childNodes[0].getElementsByTagName("a")[0]
taxid = temp.href.split("id=")[1];
taxid_name = temp.innerText;
temp.innerHTML = taxid_name + "<font color='#FF0000'>" + "(" + taxid + ")" + "</font>";
i += 1;
for (; i < cars.length; i++) {
    if (cars[i].nodeType == 1) {
        temp = cars[i].childNodes[0].getElementsByTagName("a")[0]
        taxid = temp.href.split("id=")[1];
        taxid_name = temp.innerText;
        temp.innerHTML = taxid_name + "<font color='#FF0000'>" + "(" + taxid + ")" + "</font>";
    }
}