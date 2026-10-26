'use strict';
const dataset = window.SPATIAL_DATA.samples, maxima = window.SPATIAL_DATA.vmax;
const colors = ['#3166D5','#20A4BF','#3EAD7B','#E8BB35','#E68A34','#D95577','#7754AF'];
const layers = ['L1','L2','L3','L4','L5','L6','WM'];
const $ = id => document.getElementById(id);
const SVG_NS = 'http://www.w3.org/2000/svg';
const views = [$('viewer'), $('compare-canvas')];
const geneNames = Object.keys(dataset['151673'].genes);
const zoomLevels = [1, 1.5, 2, 3, 4];
let sample = '151673', mode = 'layers', gene = 'MBP', compare = false, zoomIndex = 0, selection = null;
const viewport = new Map(views.map(c => [c, {x:0, y:0, width:1100, height:1200}]));
for (const g of geneNames) $('gene').add(new Option(g,g));
$('gene').value = gene;

function color(t) {
  const stops = [[16,14,39],[91,28,110],[170,51,126],[235,100,93],[249,212,148]];
  t = Math.min(1,Math.max(0,t));
  const n=t*4, i=Math.min(3,Math.floor(n)), r=n-i;
  return `rgb(${stops[i].map((v,k)=>Math.round(v*(1-r)+stops[i+1][k]*r)).join(',')})`;
}
function svgElement(tag, attributes={}) {
  const element=document.createElementNS(SVG_NS,tag);
  for (const [key,value] of Object.entries(attributes)) element.setAttribute(key,value);
  return element;
}
// Tissue and spots always use the same fixed coordinate transform. Zoom only
// changes the SVG viewport, never the data coordinates or expression scale.
function transform(c) {
  const width=c.id==='hero-canvas'?900:1100, height=c.id==='hero-canvas'?1000:1200, pad=45;
  const scale=Math.min((width-2*pad)/1400,(height-2*pad)/1630);
  return {scale,x:pad-330*scale,y:pad-90*scale,width,height};
}
function base(c,s,show) {
  const t=transform(c);
  c.replaceChildren(svgElement('rect',{width:t.width,height:t.height,fill:'#edf0ea'}));
  if (show) c.append(svgElement('image',{href:`assets/${s}_he.png`,x:t.x,y:t.y,width:2000*t.scale,height:2000*t.scale}));
  return t;
}
function hero() {
  const c=$('hero-canvas'),t=base(c,'151673',true),d=dataset['151673'];
  const marks=svgElement('g',{opacity:.76,'aria-hidden':'true'});
  for(let i=0;i<d.x.length;i++) marks.append(svgElement('circle',{cx:t.x+d.x[i]*t.scale,cy:t.y+d.y[i]*t.scale,r:4.7*t.scale,fill:colors[d.layer[i]]}));
  c.append(marks);
}
function viewSample(c) { return compare ? (c===views[0]?'151673':'151674') : sample; }
function highlight() {
  for(const c of views) {
    const ring=c.querySelector('.selected-spot'); if(!ring) continue;
    const active=selection && selection.view===c;
    ring.setAttribute('visibility',active?'visible':'hidden');
    if(active) {
      const t=transform(c),d=dataset[selection.sample];
      ring.setAttribute('cx',t.x+d.x[selection.index]*t.scale);
      ring.setAttribute('cy',t.y+d.y[selection.index]*t.scale);
    }
  }
}
function updateInfo() {
  $('clear-spot').hidden=!selection?.pinned;
  $('spot-state').textContent=selection?.pinned?'已固定':'移动查看';
  if(!selection) {
    $('spot-info').textContent=compare ? '151673：3,601 个捕获点\n151674：3,633 个捕获点\n点击任一切片的捕获点，固定并查看数值。' : `当前切片：${sample}\n分析点：${dataset[sample].x.length.toLocaleString()}\n点击捕获点可固定数值。`;
    return;
  }
  const d=dataset[selection.sample],i=selection.index;
  let text=`${selection.sample} · ${d.barcode[i]}\n参考层：${layers[d.layer[i]]} · UMI：${d.umi[i].toLocaleString()}\n线粒体比例：${d.pct_mt[i].toFixed(1)}%`;
  if(mode==='gene') text+=`\n${gene}：${d.genes[gene][i].toFixed(3)}`;
  if(d.pred) text+=`\n预测层：${layers[d.pred[i]]} · 置信度：${(d.confidence[i]*100).toFixed(1)}%`;
  $('spot-info').textContent=text;
}
function clearSelection() { selection=null; highlight(); updateInfo(); }
function render() {
  for(const c of views) {
    if(c===views[1] && !compare) continue;
    const s=viewSample(c),d=dataset[s],t=base(c,s,$('histology').checked||mode==='he');
    c.dataset.sample=s;
    const marks=svgElement('g',{id:c===views[0]?'spot-layer':'compare-spot-layer',opacity:Number($('opacity').value)/100,'aria-hidden':'true'});
    if(mode!=='he') for(let i=0;i<d.x.length;i++) {
      let fill;
      if(mode==='layers') fill=colors[d.layer[i]];
      else if(mode==='gene') fill=color(d.genes[gene][i]/maxima[gene]);
      else if(mode==='pred') fill=colors[d.pred[i]];
      else if(mode==='error') fill=d.pred[i]===d.layer[i]?'#b8ceca':'#d63f6b';
      else fill=color(d.entropy[i]);
      marks.append(svgElement('circle',{cx:t.x+d.x[i]*t.scale,cy:t.y+d.y[i]*t.scale,r:6.2*t.scale,fill}));
    }
    c.append(marks,svgElement('circle',{class:'selected-spot',r:10*t.scale,fill:'none',stroke:'#fff','stroke-width':2,'pointer-events':'none'}));
  }
  const left=compare?'151673':sample;
  $('slice-title').textContent=left+' · '+(left==='151673'?'开发切片':'留出切片');
  $('slice-count').textContent=dataset[left].x.length.toLocaleString()+' 个捕获点';
  const title=mode==='gene'?gene+' · log(1 + CP10K)':$('mode').selectedOptions[0].textContent;
  $('canvas-label').textContent=title; $('compare-label').textContent=title;
  highlight();
}
function update() {
  sample=$('sample').value; mode=$('mode').value; gene=$('gene').value;
  const diagnostic=['pred','error','entropy'].includes(mode);
  if(diagnostic) { compare=false; sample='151674'; $('sample').value=sample; }
  $('compare-view').disabled=diagnostic;
  $('compare-view').title=diagnostic?'预测诊断仅提供留出切片 151674':'并排显示两张切片';
  $('single-view').setAttribute('aria-pressed',String(!compare));
  $('compare-view').setAttribute('aria-pressed',String(compare));
  $('sample-control').hidden=compare;
  $('compare-panel').hidden=!compare;
  $('slice-views').classList.toggle('is-comparing',compare);
  $('gene-controls').hidden=mode!=='gene';
  $('opacity-value').textContent=$('opacity').value+'%';
  const legend=$('legend'); legend.replaceChildren();
  if(mode==='layers'||mode==='pred') {
    layers.forEach((l,i)=>{const span=document.createElement('span'),sw=document.createElement('i');sw.style.background=colors[i];span.append(sw,document.createTextNode(l));legend.append(span)});
  } else if(mode==='gene'||mode==='entropy') {
    const bar=document.createElement('div'),scale=document.createElement('div');bar.className='gradient';scale.className='scale';
    const start=document.createElement('span'),end=document.createElement('span');start.textContent='0';end.textContent=mode==='gene'?maxima[gene].toFixed(2):'1.00';scale.append(start,end);legend.append(bar,scale);
  } else if(mode==='error') legend.textContent='灰绿：正确 · 粉红：错误';
  $('mode-note').textContent=mode==='gene'?'log(1 + 每万计数)；两张切片共享色阶，上限为共同第 99 百分位。':diagnostic?'预测诊断仅提供预先留出的 151674。':mode==='he'?'H&E 组织原图；捕获点位置与该切片组织图配准。':'L1–L6 为皮层参考层，WM 为白质。仅显示质控合格且有标签的分析点。';
  selection=null; render(); updateInfo();
}
function applyViewport(c) {
  const v=viewport.get(c);
  v.x=Math.max(0,Math.min(1100-v.width,v.x));v.y=Math.max(0,Math.min(1200-v.height,v.y));
  c.setAttribute('viewBox',`${v.x} ${v.y} ${v.width} ${v.height}`);
  c.classList.toggle('is-zoomed',zoomIndex>0);
}
function zoom(next) {
  zoomIndex=Math.max(0,Math.min(zoomLevels.length-1,next));
  const z=zoomLevels[zoomIndex];
  for(const c of views) {
    const old=viewport.get(c),w=1100/z,h=1200/z;
    viewport.set(c,{x:old.x+(old.width-w)/2,y:old.y+(old.height-h)/2,width:w,height:h});applyViewport(c);
  }
  $('zoom-value').textContent=z+'×';$('zoom-out').disabled=zoomIndex===0;$('zoom-in').disabled=zoomIndex===zoomLevels.length-1;
}
function resetView() { zoom(0); }
function localPoint(c,event,matrix=c.getScreenCTM()?.inverse()) {
  if(!matrix) return null;
  const p=c.createSVGPoint();p.x=event.clientX;p.y=event.clientY;return p.matrixTransform(matrix);
}
function inspect(c,event,pin=false) {
  if(selection?.pinned && !pin) return;
  const local=localPoint(c,event); if(!local) return;
  const s=viewSample(c),d=dataset[s],t=transform(c);
  let dist=Infinity,best=-1;
  for(let i=0;i<d.x.length;i++) {const distance=(t.x+d.x[i]*t.scale-local.x)**2+(t.y+d.y[i]*t.scale-local.y)**2;if(distance<dist){dist=distance;best=i}}
  if(dist>500) { if(pin||!selection?.pinned) clearSelection(); return; }
  if(pin && selection?.pinned && selection.view===c && selection.index===best) {clearSelection();return;}
  selection={view:c,sample:s,index:best,pinned:pin};highlight();updateInfo();
}
for(const id of ['sample','mode','gene','histology']) $(id).addEventListener('change',update);
$('opacity').addEventListener('input',()=>{$('opacity-value').textContent=$('opacity').value+'%';render()});
$('single-view').addEventListener('click',()=>{compare=false;update()});
$('compare-view').addEventListener('click',()=>{compare=true;update()});
$('zoom-in').addEventListener('click',()=>zoom(zoomIndex+1));
$('zoom-out').addEventListener('click',()=>zoom(zoomIndex-1));
$('reset-view').addEventListener('click',resetView);
$('clear-spot').addEventListener('click',clearSelection);
for(const c of views) {
  let drag=null;
  c.addEventListener('pointerdown',event=>{
    if(event.button!==0 || event.isPrimary===false) return;
    drag={id:event.pointerId,x:event.clientX,y:event.clientY,view:{...viewport.get(c)},matrix:c.getScreenCTM()?.inverse(),moved:false};
    if(zoomIndex>0) c.setPointerCapture(event.pointerId);
  });
  c.addEventListener('pointermove',event=>{
    if(drag && event.buttons===0) {drag=null;c.classList.remove('is-dragging')}
    if(drag && event.pointerId===drag.id) {
      if(Math.hypot(event.clientX-drag.x,event.clientY-drag.y)>5) drag.moved=true;
      if(drag.moved && zoomIndex>0) {
        const current=localPoint(c,event,drag.matrix),start=localPoint(c,{clientX:drag.x,clientY:drag.y},drag.matrix);
        if(current&&start) {viewport.set(c,{...drag.view,x:drag.view.x-(current.x-start.x),y:drag.view.y-(current.y-start.y)});applyViewport(c);c.classList.add('is-dragging')}
        return;
      }
    }
    if(event.pointerType!=='touch') inspect(c,event);
  });
  c.addEventListener('pointerup',event=>{
    if(!drag || event.pointerId!==drag.id) return;
    const click=!drag.moved;drag=null;c.classList.remove('is-dragging');
    if(c.hasPointerCapture(event.pointerId))c.releasePointerCapture(event.pointerId);
    if(click)inspect(c,event,true);
  });
  c.addEventListener('pointercancel',()=>{drag=null;c.classList.remove('is-dragging')});
  c.addEventListener('pointerleave',()=>{if(!drag&&!selection?.pinned)clearSelection()});
  c.addEventListener('keydown',event=>{
    if(event.key==='+'||event.key==='=')zoom(zoomIndex+1);
    else if(event.key==='-')zoom(zoomIndex-1);
    else if(event.key==='0')resetView();
    else if(event.key==='Escape')clearSelection();
    else if(['ArrowLeft','ArrowRight','ArrowUp','ArrowDown'].includes(event.key)&&zoomIndex>0) {
      const v=viewport.get(c);v.x+=event.key==='ArrowLeft'?-v.width*.12:event.key==='ArrowRight'?v.width*.12:0;v.y+=event.key==='ArrowUp'?-v.height*.12:event.key==='ArrowDown'?v.height*.12:0;applyViewport(c);
    } else return;
    event.preventDefault();
  });
}
zoom(0);

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
