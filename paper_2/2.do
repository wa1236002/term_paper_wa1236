clear all
cls
set more off


* 数据导入与处理
import excel "C:\Users\lenovo\Desktop\area_price.xlsx", sheet("Sheet1") firstrow clear
rename 编号 id
rename 总面积 area
rename 价格 price
drop if missing(area) | missing(price)
summarize id area price

* 核密度估计
* 核密度估计图
kdensity price, kernel(gauss) normal ///
    title("上海二手房价格高斯核密度估计（最优带宽）") ///
    xtitle("价格（万元）") ytitle("密度 (Density)") ///
    name(kde_full, replace)
* 0-3000万主体区间图
kdensity price if price <= 3000, kernel(gauss) normal ///
    title("上海二手房价格高斯核密度估计 (最优带宽,局部)") ///
    xtitle("价格（万元）") ytitle("密度 (Density)") ///
    name(kde_local, replace)
* 带宽敏感性对比实验：欠平滑h=10
kdensity price if price <= 3000, kernel(gauss) normal bwidth(10) ///
    title("上海二手房价格高斯核密度估计 (h=10欠平滑)") ///
    xtitle("价格（万元）") ytitle("密度 (Density)") ///
    name(kde_h10, replace)
* 带宽敏感性对比实验：过平滑h=50
kdensity price if price <= 3000, kernel(gauss) normal bwidth(50) ///
    title("上海二手房价格高斯核密度估计 (h=50过平滑)") ///
    xtitle("价格（万元）") ytitle("密度 (Density)") ///
    name(kde_h50, replace)


* 非参数回归分析
* 非参数回归
lpoly price area, noscatter kernel(gauss) ///
    title("全样本非参数回归") ///
    xtitle("总面积 (平方米)") ytitle("价格 (万元)") ///
    name(lpoly_auto, replace)
* 扩大压轴带宽（h=30）
lpoly price area, noscatter kernel(gauss) bwidth(30) ///
    title("全样本非参数回归 (大带宽h=30)") ///
    xtitle("总面积 (平方米)") ytitle("价格 (万元)") ///
    name(lpoly_30, replace)

	
* Bootstrap重抽样与置信区间估计
lpoly price area, noscatter kernel(gauss) ci reps(500) seed(1236) ///
    title("全样本非参数回归及Bootstrap置信区间") ///
    xtitle("总面积 (平方米)") ytitle("价格 (万元)") ///
    name(lpoly_bootstrap_auto, replace)
	
	
* Wilcoxon秩和检验
gen size_group = (area > 100)
twoway (kdensity price if size_group==0, kernel(gauss) color(blue) lpattern(solid)) ///
       (kdensity price if size_group==1, kernel(gauss) color(red) lpattern(dash)) ///
       if price <= 3000, ///
       title("基于2024新政标准的刚需房与改善房价格KDE对比") ///
       legend(label(1 "新政刚需型 (<=100㎡)") label(2 "新政改善型 (>100㎡)")) ///
       xtitle("价格 (万元)") ytitle("密度") name(kde_compare_100, replace)
* 运行 Wilcoxon 秩和检验
ranksum price, by(size_group)
