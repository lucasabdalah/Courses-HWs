<!DOCTYPE html>
<html>
<script>
function includeHTML() {
  var z, i, elmnt, file, xhttp;
  /*loop through a collection of all HTML elements:*/
  z = document.getElementsByTagName("*");
  for (i = 0; i < z.length; i++) {
    elmnt = z[i];
    /*search for elements with a certain atrribute:*/
    file = elmnt.getAttribute("w3-include-html");
    if (file) {
      /*make an HTTP request using the attribute value as the file name:*/
      xhttp = new XMLHttpRequest();
      xhttp.onreadystatechange = function() {
        if (this.readyState == 4) {
          if (this.status == 200) {elmnt.innerHTML = this.responseText;}
          if (this.status == 404) {elmnt.innerHTML = "Page not found.";}
          /*remove the attribute, and call this function once more:*/
          elmnt.removeAttribute("w3-include-html");
          includeHTML();
        }
      }      
      xhttp.open("GET", file, true);
      xhttp.send();
      /*exit the function:*/
      return;
    }
  }
};
</script>
<h1>
    <center> 
        Lucas Abdalah <br>
        [TI8419 - Multilinear Algebra] Homeworks <br>
        Professors: André Lima e Henrique Goulart <br>
    </center>
</h1>
<div id="toc">
    <h1> Contents </h1>
    <a href="#hw0-report">HW0 report </a> <br>
    <a href="#hw1-report">HW1 report </a> <br>
    <a href="#hw2-report">HW2 report </a> <br>
    <a href="#hw3-report">HW3 report </a> <br>
    <a href="#hw4-report">HW4 report </a> <br>
    <a href="#hw5-report">HW5 report </a> <br>
    <a href="#hw6-report">HW6 report </a> <br>
    <a href="#hw7-report">HW7 report </a> <br>
    <a href="#hw8-report">HW8 report </a> <br>
    <a href="#hw9-report">HW9 report </a> <br>
    <a href="#hw10-report">HW10 report </a> <br>
    <a href="#hw11-report">HW11 report </a> <br>
    <a href="#hw12-report">HW12 report </a> <br>
    <a href="#hw13-report">HW13 report </a> <br>
    <a href="#code-report">Code report </a> <br>
</div>
<hr>
<br>
<div> 
<div id="hw0-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw0/hw0-report.html"></div>
<br> 
<div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> Return to Table of Contents </a> 
    </p> 
</div> 
<br>
<div id="hw1-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw1/hw1-report.html">
</div>
<br> 
<div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div> 
<br>
<div id="hw2-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw2/hw2-report.html"></div>
<br> 
<div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div> 
<br>
<div id="hw3-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw3/hw3-report.html"></div>
<br> 
<div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div> 
<br>
<div id="hw4-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw4/hw4-report.html"></div>
<br> 
<div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div> 
<br>
<div id="hw5-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw5/hw5-report.html"></div>
<br> 
<div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div> 
<br>
<div id="hw6-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw6/hw6-report.html"></div>
<br> 
<div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div> 
<br>
<div id="hw7-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw7/hw7-report.html"></div>
<br> <div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div> 
<br>
<div id="hw8-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw8/hw8-report.html"></div>
<br> <div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div> 
<br>
<div id="hw9-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw9/hw9-report.html"></div>
<br> <div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div>
<br>
<div id="hw10-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw10/hw10-report.html"></div>
<br> <div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div>
<br>
<div id="hw11-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw11/hw11-report.html"></div>
<br> <div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div>
<br>
<div id="hw12-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw12/hw12-report.html"></div>
<br> <div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div>
<br>
<div id="hw13-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/hw13/hw13-report.html"></div>
<br> <div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div>
<div id="code-report" w3-include-html="https://raw.githubusercontent.com/lucasabdalah/Courses-HWs/algebra/Master/TIP8419-ALGEBRA_LINEAR_E_MULTILINEAR/homework/nd-web.html"></div>
<br> <div style="background-color:rgba(200, 0, 000, 0.15); text-align:justify; padding:1px"> 
    <p> 
        <a href="#toc"> 
            Return to Table of Contents 
        </a> 
    </p> 
</div>
<br>
</div>
<script>
includeHTML();
</script>
</body>
</html>
