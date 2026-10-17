'use strict';
const dataset=window.SPATIAL_DATA.samples, maxima=window.SPATIAL_DATA.vmax;
const colors=['#3166D5','#20A4BF','#3EAD7B','#E8BB35','#E68A34','#D95577','#7754AF'];
const layers=['L1','L2','L3','L4','L5','L6','WM'];
const $=id=>document.getElementById(id), canvas=$('viewer');
const SVG_NS='http://www.w3.org/2000/svg';
const geneNames=Object.keys(dataset['151673'].genes);
let sample='151673',mode='layers',gene='MBP', selected=-1;
for(const g of geneNames){const o=new Option(g,g);$('gene').add(o)}$('gene').value=gene;
function color(t){const stops=[[16,14,39],[91,28,110],[170,51,126],[235,100,93],[249,212,148]];t=Math.min(1,Math.max(0,t));const n=t*4,i=Math.min(3,Math.floor(n)),r=n-i;return `rgb(${stops[i].map((v,k)=>Math.round(v*(1-r)+stops[i+1][k]*r)).join(',')})`}
// Vector marks share the same image-pixel mapping as the stored coordinates.
function svgElement(tag,attributes={}){
  const element=document.createElementNS(SVG_NS,tag);
  for(const [key,value] of Object.entries(attributes))element.setAttribute(key,value);
  return element;
}
function transform(c){
  const {width,height}=c.viewBox.baseVal,pad=45;
  const scale=Math.min((width-2*pad)/1400,(height-2*pad)/1630);
  return {scale,x:pad-330*scale,y:pad-90*scale};
}
function base(c,s,show){
  const t=transform(c),{width,height}=c.viewBox.baseVal;
  c.replaceChildren(svgElement('rect',{width,height,fill:'#edf0ea'}));
  if(show)c.append(svgElement('image',{href:`assets/${s}_he.png`,x:t.x,y:t.y,width:2000*t.scale,height:2000*t.scale}));
  return t;
}
function hero(){
  const c=$('hero-canvas'),t=base(c,'151673',true),d=dataset['151673'];
  const marks=svgElement('g',{opacity:.76,'aria-hidden':'true'});
  for(let i=0;i<d.x.length;i++)marks.append(svgElement('circle',{cx:t.x+d.x[i]*t.scale,cy:t.y+d.y[i]*t.scale,r:4.7*t.scale,fill:colors[d.layer[i]]}));
  c.append(marks);
}
function highlight(){
  const ring=$('selected-spot');if(!ring)return;
  ring.setAttribute('visibility',selected<0?'hidden':'visible');
  if(selected>=0){const t=transform(canvas),d=dataset[sample];ring.setAttribute('cx',t.x+d.x[selected]*t.scale);ring.setAttribute('cy',t.y+d.y[selected]*t.scale);}
}
function render(){
  const d=dataset[sample],t=base(canvas,sample,$('histology').checked||mode==='he');
  const marks=svgElement('g',{id:'spot-layer',opacity:Number($('opacity').value)/100,'aria-hidden':'true'});
  if(mode!=='he')for(let i=0;i<d.x.length;i++){
    let fill;
    if(mode==='layers')fill=colors[d.layer[i]];
    else if(mode==='gene')fill=color(d.genes[gene][i]/maxima[gene]);
    else if(mode==='pred')fill=colors[d.pred[i]];
    else if(mode==='error')fill=d.pred[i]===d.layer[i]?'#b8ceca':'#d63f6b';
    else fill=color(d.entropy[i]);
    marks.append(svgElement('circle',{cx:t.x+d.x[i]*t.scale,cy:t.y+d.y[i]*t.scale,r:6.2*t.scale,fill}));
  }
  canvas.append(marks,svgElement('circle',{id:'selected-spot',r:10*t.scale,fill:'none',stroke:'#fff','stroke-width':2,'pointer-events':'none'}));
  highlight();
  $('canvas-label').textContent=sample+' · '+(mode==='gene'?gene+' 表达':$('mode').selectedOptions[0].textContent);
}
function update(){sample=$('sample').value;mode=$('mode').value;gene=$('gene').value;selected=-1;const diagnostic=['pred','error','entropy'].includes(mode);if(diagnostic&&sample==='151673'){sample='151674';$('sample').value=sample}$('gene-controls').hidden=mode!=='gene';$('opacity-value').textContent=$('opacity').value+'%';const legend=$('legend');legend.replaceChildren();if(mode==='layers'||mode==='pred'){layers.forEach((l,i)=>{const span=document.createElement('span');const sw=document.createElement('i');sw.style.background=colors[i];span.append(sw,document.createTextNode(l));legend.append(span)})}else if(mode==='gene'||mode==='entropy'){const bar=document.createElement('div');bar.className='gradient';const scale=document.createElement('div');scale.className='scale';const end=mode==='gene'?maxima[gene].toFixed(2):'1.00';scale.innerHTML=`<span>0</span><span>${end}</span>`;legend.append(bar,scale)}else if(mode==='error'){legend.textContent='灰绿：正确 · 粉红：错误'}$('mode-note').textContent=mode==='gene'?'log(1 + 每万计数)；同一基因在两张切片使用相同色阶，截断于共同第 99 百分位。':diagnostic?'预测诊断只对预先留出的 151674 提供，选择该模式时自动切换。':'L1–L6 为皮层参考层，WM 为白质。仅显示质控合格且有标签的分析点。';$('spot-info').textContent=`当前切片：${sample}\n分析点：${dataset[sample].x.length.toLocaleString()}\n移动或触摸一个捕获点查看数值。`;render()}
for(const id of ['sample','mode','gene','histology'])$(id).addEventListener('change',update);$('opacity').addEventListener('input',()=>{$('opacity-value').textContent=$('opacity').value+'%';render()});
function inspect(event){
  const point=canvas.createSVGPoint();point.x=event.clientX;point.y=event.clientY;
  const matrix=canvas.getScreenCTM();if(!matrix)return;
  const local=point.matrixTransform(matrix.inverse()),t=transform(canvas),d=dataset[sample];
  let dist=Infinity,best=-1;
  for(let i=0;i<d.x.length;i++){const z=(t.x+d.x[i]*t.scale-local.x)**2+(t.y+d.y[i]*t.scale-local.y)**2;if(z<dist){dist=z;best=i;}}
  if(dist>500){selected=-1;highlight();return;}
  selected=best;
  let text=`${d.barcode[best]}\n参考层：${layers[d.layer[best]]} · UMI：${d.umi[best]}\n线粒体比例：${d.pct_mt[best].toFixed(1)}%`;
  if(mode==='gene')text+=`\n${gene}：${d.genes[gene][best].toFixed(3)}`;
  if(d.pred)text+=`\n预测层：${layers[d.pred[best]]} · 置信度：${(d.confidence[best]*100).toFixed(1)}%`;
  $('spot-info').textContent=text;highlight();
}
canvas.addEventListener('pointermove',inspect);canvas.addEventListener('pointerdown',inspect);
const names={Majority:'多数类基线',Logistic:'逻辑回归',RandomForest:'随机森林'};
for(const r of window.MODEL_RESULTS){const tr=document.createElement('tr');for(const v of [names[r.model],r.accuracy.toFixed(3),r.macro_f1.toFixed(3),r.macro_auroc.toFixed(3),r.macro_ap.toFixed(3)]){const td=document.createElement('td');td.textContent=v;tr.append(td)}$('metrics').append(tr)}
const figs=[['01_mindmap','分析流程与验证设计','开发切片拟合，留出切片评估'],['02_histology_layers','原始切片与参考分层','两张 H&E 与配准叠加'],['03_quality_control','质量控制','深度、基因数与线粒体比例'],['04_layer_composition','各层组成与质量分布','样本不均衡与分布差异'],['05_marker_spatial_maps','标志基因空间表达','六基因、两切片共同色阶'],['06_heatmap_moran','表达热图与空间自相关','层均值和 Moran I'],['07_layer_marker_boxplots','层间表达分布','CUX2 / PCP4 / MBP'],['08_embeddings','PCA 与 t-SNE','表达结构与参考层'],['09_clustering_comparison','聚类方法对照','表达聚类与邻域平滑'],['10_model_comparison','模型与验证设计','随机、空间块与跨切片'],['11_confusion_matrices','混淆矩阵','各层的预测分配'],['12_roc_pr','ROC 与 PR 曲线','排序能力与类别不平衡'],['13_calibration_uncertainty','校准与不确定性','最大概率、正确率与预测熵'],['14_spatial_prediction_errors','空间误差定位','参考层、预测与错误分布'],['15_interpretation','模型解释与统计关联','回投系数与稳健区间']];
figs.push(...[["17_histology_multiscale", "组织全貌与局部形态", "原始 H&E 定位裁剪与成对视野"], ["18_cross_section_replication", "跨切片表达一致性", "2,000 个基因的层对比"], ["19_marker_discovery_validation", "候选基因发现与复核", "开发片筛选，留出片检查方向"], ["20_spatial_expression_gradients", "连续空间表达梯度", "几何深度与 CUX2 / RORB / MBP"], ["21_layer_neighborhood_network", "参考层邻域接触", "置换富集与层间接触网络"], ["22_multiscale_spatial_structure", "多尺度空间结构", "12 个基因 × 5 种邻域半径"], ["23_clustering_stability", "聚类稳定性", "5 种权重与 8 个随机种子"], ["24_prediction_error_structure", "错误层级与边界梯度", "选择性预测和空间块区间"], ["25_model_disagreement", "两种模型的分歧", "预测分歧的空间分布"], ["26_image_expression_ablation", "图像与表达消融", "三种输入，同一留出切片"]]);
for(const [file,title,sub] of figs){const button=document.createElement('button');button.type='button';const img=document.createElement('img');img.src=`assets/${file}.svg`;img.alt=title;img.loading='lazy';const b=document.createElement('b');b.textContent=title;const span=document.createElement('span');span.textContent=sub;button.append(img,b,span);button.addEventListener('click',()=>{$('dialog-image').src=img.src;$('dialog-image').alt=title;$('dialog-caption').textContent=title+' · '+sub;$('download-svg').href=`assets/${file}.svg`;$('image-dialog').showModal()});$('figure-grid').append(button)}
$('close-dialog').addEventListener('click',()=>$('image-dialog').close());$('image-dialog').addEventListener('click',e=>{if(e.target===$('image-dialog'))$('image-dialog').close()});
for(const r of window.REFERENCES){const li=document.createElement('li');const a=document.createElement('a');a.href='https://doi.org/'+r.doi;a.target='_blank';a.rel='noopener';a.textContent=`${r.title}. ${r.journal}, ${r.year}.`;li.append(a);$('references').append(li)}
hero();
update();

const modalNames={'Image only':'仅图像（8 维）',Logistic:'仅表达（50 维）','Expression + image':'表达 + 图像（58 维）'};
for(const name of Object.keys(modalNames)){const r=window.EXTENDED_RESULTS.find(x=>x.model===name),tr=document.createElement('tr');for(const v of [modalNames[name],r.accuracy.toFixed(3),r.macro_f1.toFixed(3),r.macro_auroc.toFixed(3)]){const td=document.createElement('td');td.textContent=v;tr.append(td)}$('modality-metrics').append(tr)}
