(function(mod) {
  if (typeof exports == "object" && typeof module == "object") // CommonJS
    mod(require("../../lib/codemirror"), require("../../addon/mode/simple"));
  else if (typeof define == "function" && define.amd) // AMD
    define(["../../lib/codemirror", "../../addon/mode/simple"], mod);
  else // Plain browser env
    mod(CodeMirror);
})(function(CodeMirror) {
"use strict";


CodeMirror.defineSimpleMode("dacalc", {
    // The start state contains the rules that are intially used
    start: [
	{regex: /\#[^\n]*/, token: "comment"},
	{regex: /\[[^\]]*\]/, token: "string"},
	{regex: /[\'\"][^\'\"]*[\'\"]/, token: "builtin"},// maybe find better token name some time
	{regex: /def|use|output|dim|import|analyze|\?/, token: "keyword"},
	{regex: /sqrt|sin|cos|tan|asin|acos|atan|atan2|abs|log|ln|log2|log10|pow/, token: "keyword"},
	{regex: /_[a-zA-Z0-9_]+/, token: "atom"},
        {regex: /[a-zA-Z][a-zA-Z0-9_]*/, token: "variable"},
	{regex: /(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?/, token: "number"},
	{regex: /[\+\-\*/=\^\%]+/, token: "operator"}
    ],
  meta: {
    lineComment: "\#"
  }
});

});

