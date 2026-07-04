* 删除已加载的数据并清屏
clear
cls
* 数据准备与重命名
graph set window fontface "SimSun"
import excel "C:\Users\lenovo\Desktop\stata_ready.xlsx", sheet("Sheet1") firstrow
describe
encode 地区, gen(province_id)
xtset province_id 年份
rename 全体居民人均可支配收入元 x1
rename 全体居民人均消费支出元 x2
rename 每万人拥有医疗卫生机构数个万人 x3
rename 每万人拥有医疗卫生床位数张万人 x4
rename 每万人拥有卫生技术人员数人 x5
rename 每万人拥有执业助理医师数人 x6
rename 千人拥有本专科在校生数人千人 x7
rename 千人拥有高中在校生数人千人 x8
rename 每万人专利申请授权量项万人 x9
rename 每万人发明专利申请授权量项万人 x10
rename 平均每人拥有移动电话量部人 x11
rename 人均快递包裹数件人 x12
rename 城市绿地空间分布密度 x13
rename 建成区绿化覆盖率 x14
rename 路网密度万平方米平方公里 x15
rename 城市公厕分布密度座平方公里 x16
rename 城市道路保洁覆盖密度 x17

* 描述性统计
summarize x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17

* 主成分分析
pca x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17

* 因子分析
factortest x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17
factor x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17, pf factors(4)
rotate, varimax
capture drop f1 f2 f3 f4
predict f1 f2 f3 f4
capture drop comp_score
gen comp_score = (5.73148*f1 + 2.93111*f2 + 1.51357*f3 + 1.34973*f4) / 11.52589
label variable comp_score "城市综合发展水平得分"
capture drop rank
bysort 年份: egen rank = rank(comp_score), field
sort 年份 rank
browse 地区 年份 f1 f2 f3 f4 comp_score rank

* 聚类分析
* 聚类谱系图
preserve
keep if 年份 == 2024
capture cluster drop my_cluster
cluster wardslinkage f1 f2 f3 f4, name(my_cluster)
cluster dendrogram my_cluster, ///
    label(地区) ///
    xlabel(, angle(90) labsize(small)) ///
    title("2024年各省份城市综合发展水平聚类谱系图") ///
    ytitle("平方欧氏距离 (Dissimilarity)") ///
    name(g1, replace)
capture drop group_id
cluster generate group_id = groups(5), name(my_cluster)
sort group_id comp_score
list 地区 comp_score group_id
* 划分散点图
capture drop new_group
egen new_group = group(group_id)
capture drop comp_score1 comp_score2 comp_score3 comp_score4 comp_score5
separate comp_score, by(new_group)
twoway (scatter comp_score1 rank, mlabel(地区) msize(small) mlabangle(45) mlabsize(vsmall)) ///
       (scatter comp_score2 rank, mlabel(地区) msize(small) mlabangle(45) mlabsize(vsmall)) ///
       (scatter comp_score3 rank, mlabel(地区) msize(small) mlabangle(45) mlabsize(vsmall)) ///
       (scatter comp_score4 rank, mlabel(地区) msize(small) mlabangle(45) mlabsize(vsmall)) ///
       (scatter comp_score5 rank, mlabel(地区) msize(small) mlabangle(45) mlabsize(vsmall)), ///
       title("2024年各省份城市综合发展水平梯队划分与排名分布") ///
       xtitle("全国综合发展水平排名 (Rank)") ///
       ytitle("城市综合发展水平得分 (Score)") ///
       legend(label(1 "第一梯队") label(2 "第二梯队") label(3 "第三梯队") label(4 "第四梯队") label(5 "第五梯队")) ///
       xline(1(1)31, lstyle(grid) lcolor(gainsboro*0.5)) ///
       yline(-1(0.5)3, lstyle(grid) lcolor(gainsboro*0.5)) ///
       name(g2, replace)

graph display g1
graph display g2
restore
