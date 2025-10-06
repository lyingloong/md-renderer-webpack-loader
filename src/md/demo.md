# 测试 Markdown#2025-03-05
#fig
这是一张图片：
![a figure](https://static.igem.wiki/teams/5924/drylab/snap-25.webp)

# basic inline elements and math formulas
some *text* with $\LaTeX$ ^here^, which _may be_ very *^very^* long _^long^_, *even wraps-the-line*
\[ display\ mode\ math\ and\ huge\ operators\ like \sum_{i=1}^n \]
and another line
#link and figure
[and a link](/)
#lists
+ test
	Plain text
	+ another element
		+ and embedded element
			embedded text
		+ aligned element
+ txt

Outer

#display math wrapped and $inline~\LaTeX~in~title$
\[
	display\ mode\ math\ that\ wraps\ the\ line\ in\ source\ code
\]
#code-block
	```javascript=2
	console.log('hello,\n world')
	```
Outer
#table
| Name | Age | City |
| ---- | ---:| :---:|
| Tom  |  12 | NY   |
| Lily |  22 | LA   |