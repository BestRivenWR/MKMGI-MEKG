/* Real source-backed relation explorer; project records stay in this browser. */
(() => {
 'use strict';
 const ns='http://www.w3.org/2000/svg';
 const make=(tag,attrs={},label)=>{const n=document.createElementNS(ns,tag);for(const [k,v] of Object.entries(attrs))n.setAttribute(k,String(v));if(label!==undefined)n.textContent=label;return n;};
 const colors={classic_experimental:'#b58228',literature_reported:'#b99564',public_experimental:'#95752f',public_association:'#8775a2',ppi:'#678797',model_regulatory:'#519c83',project_inferred:'#1d7659'};
 const angles={classic_experimental:-90,literature_reported:-40,public_experimental:10,ppi:60,public_association:110,project_inferred:160,model_regulatory:215};
 class Explorer {
  constructor(host,options){this.host=host;this.options=options;this.nodes=new Map();this.edges=new Map();this.layers=new Set(Object.keys(colors));this.root=null;this.transform={x:0,y:0,k:1};}
  reset(gene){this.root=gene;this.nodes.clear();this.edges.clear();this.layers=new Set(Object.keys(colors));this.host.replaceChildren();}
  add(data){data.nodes.forEach(n=>this.nodes.set(n.id,n));data.edges.forEach(e=>this.edges.set(e.id,e));this.draw();}
  filter(layer,yes){if(yes)this.layers.add(layer);else this.layers.delete(layer);this.draw();}
  zoom(by){const old=this.transform.k;this.transform.k=Math.max(.04,Math.min(8,old*by));this.transform.x=520-(520-this.transform.x)*this.transform.k/old;this.transform.y=260-(260-this.transform.y)*this.transform.k/old;this.apply();}
  apply(){if(this.world)this.world.setAttribute('transform',`translate(${this.transform.x} ${this.transform.y}) scale(${this.transform.k})`);}
  fit(){if(!this.bounds)return;const b=this.bounds,k=Math.min(1.1,960/(b.maxx-b.minx+130),450/(b.maxy-b.miny+100));this.transform={k,x:520-k*(b.minx+b.maxx)/2,y:260-k*(b.miny+b.maxy)/2};this.apply();}
  draw(){
   const all=[...this.edges.values()].filter(e=>this.layers.has(e.layer));const active=new Set([this.root]);all.forEach(e=>{active.add(e.source);active.add(e.target);});
   const label=(en,zh)=>window.MEKGLabel?window.MEKGLabel(en,zh):en;
   this.host.replaceChildren();if(all.length===0){const p=document.createElement('p');p.className='muted';p.textContent=label('No relationships in this view.','当前视图暂无关系。');this.host.append(p);return;}
   const svg=make('svg',{viewBox:'0 0 1040 520',role:'img','aria-label':label('Gene network: published evidence, PPI and model-input relationships','基因网络：已发表证据、PPI 与模型输入关系')});this.svg=svg;
   const defs=make('defs');for(const [layer,c] of Object.entries(colors)){const marker=make('marker',{id:'arrow-'+layer,markerWidth:7,markerHeight:7,refX:7,refY:3.5,orient:'auto',markerUnits:'userSpaceOnUse'});marker.append(make('path',{d:'M0 0L7 3.5L0 7Z',fill:c}));defs.append(marker);}svg.append(defs);
   const world=make('g');this.world=world;svg.append(world);const pos=new Map([[this.root,{x:0,y:0}]]),groups={};
   const votes=new Map();for(const e of all){for(const id of [e.source,e.target]){if(id===this.root)continue;if(!votes.has(id))votes.set(id,{});const v=votes.get(id);v[e.layer]=(v[e.layer]||0)+1;}}
   for(const id of active){if(id===this.root)continue;const v=votes.get(id)||{};const dominant=Object.keys(v).sort((a,b)=>v[b]-v[a])[0]||'ppi';(groups[dominant]??=[]).push(id);}
   for(const [layer,ids] of Object.entries(groups)){
    ids.sort((a,b)=>(this.nodes.get(a)?.label||a).localeCompare(this.nodes.get(b)?.label||b));
    let index=0,ring=0;while(index<ids.length){const capacity=4+ring*2,chunk=ids.slice(index,index+capacity),r=170+ring*100;
     chunk.forEach((id,j)=>{const theta=(angles[layer]-20+(chunk.length===1?20:40*j/(chunk.length-1)))*Math.PI/180;pos.set(id,{x:r*Math.cos(theta),y:r*Math.sin(theta),layer});});index+=chunk.length;ring++;}
   }
   const xs=[...pos.values()].map(p=>p.x),ys=[...pos.values()].map(p=>p.y);this.bounds={minx:Math.min(...xs),maxx:Math.max(...xs),miny:Math.min(...ys),maxy:Math.max(...ys)};
   const connected=new Map();for(const e of all){const a=pos.get(e.source),b=pos.get(e.target);if(!a||!b)continue;
    const sign=e.id.charCodeAt(e.id.length-1)%2?1:-1,dx=b.x-a.x,dy=b.y-a.y,len=Math.hypot(dx,dy)||1,curve=sign*9;
    const startRadius=e.source===this.root?37:15,endRadius=e.target===this.root?38:16;
    const path=make('path',{d:`M${a.x+dx/len*startRadius} ${a.y+dy/len*startRadius}Q${(a.x+b.x)/2-dy/len*curve} ${(a.y+b.y)/2+dx/len*curve} ${b.x-dx/len*endRadius} ${b.y-dy/len*endRadius}`,fill:'none',stroke:colors[e.layer],'stroke-width':['classic_experimental','public_experimental'].includes(e.layer)?2.3:1.35,opacity:['classic_experimental','public_experimental'].includes(e.layer)?.9:.55,class:'network-edge',tabindex:0,role:'button','aria-label':e.source+' '+e.predicate+' '+e.target});
    if(['literature_reported','public_association','model_regulatory','project_inferred'].includes(e.layer))path.setAttribute('stroke-dasharray',e.layer==='project_inferred'?'6 3':'3 4');
    if(e.directed)path.setAttribute('marker-end','url(#arrow-'+e.layer+')');
    const edgeGroup=make('g',{tabindex:0,role:'button','aria-label':`${e.source} ${e.predicate} ${e.target}`,class:'network-edge-target'});
    for(const attr of ['tabindex','role','aria-label'])path.removeAttribute(attr);
    const hit=make('path',{d:path.getAttribute('d'),fill:'none',stroke:'transparent','stroke-width':16,'pointer-events':'stroke','aria-hidden':'true'});
    edgeGroup.append(make('title',{},`${e.source} ${e.predicate} ${e.target}`),hit,path);
    edgeGroup.addEventListener('click',ev=>{ev.stopPropagation();this.selectEdge(e,path);});edgeGroup.addEventListener('keydown',ev=>{if(ev.key==='Enter'||ev.key===' '){ev.preventDefault();this.selectEdge(e,path);}});world.append(edgeGroup);
    for(const id of [e.source,e.target]){if(!connected.has(id))connected.set(id,[]);connected.get(id).push(path);}
   }
   for(const [id,p] of pos){const n=this.nodes.get(id)||{id,label:id};const focus=id===this.root,r=focus?35:13,g=make('g',{transform:`translate(${p.x} ${p.y})`,tabindex:0,role:'button',class:'network-node','aria-label':n.label+' '+id});
    if(focus)g.append(make('circle',{r:46,fill:'#e8f1eb',stroke:'#d2e4d7','stroke-width':1}));
    g.append(make('circle',{r,fill:focus?'#196b50':'#fff',stroke:focus?'#196b50':colors[p.layer],'stroke-width':focus?2:2.4}));
    if(!focus)g.append(make('circle',{r:4,fill:colors[p.layer],opacity:.8}));
    const display=n.label===id?'v5:'+id.replace('Zm00001eb',''):String(n.label).toLowerCase()==='o2'?'O2':n.label;
    g.append(make('text',{x:0,y:focus?5:29,'text-anchor':'middle','font-size':focus?16:10,'font-weight':focus?650:500,fill:focus?'#fff':'#344e44'},display.length>18?display.slice(0,16)+'…':display));
    if(focus)g.append(make('text',{x:0,y:64,'text-anchor':'middle','font-size':11,fill:'#196b50'},label('Query gene','查询基因')));
    g.append(make('title',{},`${n.label}\n${id}\n${n.description||''}`));
    const choose=()=>{this.options.onNode(n,(connected.get(id)||[]).length);};g.addEventListener('click',ev=>{ev.stopPropagation();choose();});g.addEventListener('keydown',ev=>{if(ev.key==='Enter')choose();});
    g.addEventListener('pointerenter',()=>{world.querySelectorAll('.network-edge').forEach(x=>x.style.opacity='.12');(connected.get(id)||[]).forEach(x=>x.style.opacity='.95');});
    g.addEventListener('pointerleave',()=>world.querySelectorAll('.network-edge').forEach(x=>x.style.opacity=''));
    world.append(g);
   }
   this.host.append(svg);this.fit();let drag=null;
   svg.addEventListener('pointerdown',ev=>{if(ev.target.closest('.network-node')||ev.target.closest('.network-edge-target'))return;drag={x:ev.clientX,y:ev.clientY,tx:this.transform.x,ty:this.transform.y};svg.setPointerCapture(ev.pointerId);});
   svg.addEventListener('pointermove',ev=>{if(!drag)return;const scale=1040/svg.getBoundingClientRect().width;this.transform.x=drag.tx+(ev.clientX-drag.x)*scale;this.transform.y=drag.ty+(ev.clientY-drag.y)*scale;this.apply();});
   svg.addEventListener('pointerup',()=>{drag=null;});svg.addEventListener('wheel',ev=>{if(!ev.ctrlKey&&!ev.shiftKey)return;ev.preventDefault();this.zoom(ev.deltaY<0?1.15:1/1.15);},{passive:false});
  }
  selectEdge(e,path){this.world.querySelectorAll('.network-edge').forEach(x=>x.classList.remove('selected-edge'));path.classList.add('selected-edge');this.options.onEdge(e,this.nodes);}
 }
 window.MEKGNetwork=Explorer;
})();
