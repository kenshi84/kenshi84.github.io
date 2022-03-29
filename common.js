function switch_language(language) {
    var elements;
    elements = document.getElementsByClassName("ja"); for (var i = 0; i < elements.length; ++i) elements[i].style.display = language == "ja" ? "" : "none";
    elements = document.getElementsByClassName("en"); for (var i = 0; i < elements.length; ++i) elements[i].style.display = language == "en" ? "" : "none";
}
function init() {
    var header = document.getElementById("header");
    if (header != null)
        header.innerHTML =
            "<button class='en' onclick='switch_language(\"ja\")'>日本語 </button>" +
            "<button class='ja' onclick='switch_language(\"en\")'>English</button>";
    var footer = document.getElementById("footer");
    if (footer != null)
        footer.innerHTML =
            "<hr />" +
            "<div>&copy;<a href='http://research.nii.ac.jp/~takayama/'>Kenshi Takayama</a></div>" +
            "<div>Last modified: " + document.lastModified + "</div>";
    var language = (window.navigator.userLanguage || window.navigator.language || window.navigator.browserLanguage).substr(0,2) == "ja" ? "ja" : "en";
    switch_language(language);
}
