"""Read native VTK output independently and exercise failure/size guards."""
import argparse,json,os,re,struct,subprocess,tempfile
from pathlib import Path
import xml.etree.ElementTree as ET


def xml_arrays(path):
    data=path.read_bytes();end=data.index(b'<AppendedData')
    header=ET.fromstring(data[:end]+b'</VTKFile>')
    start=data.index(b'_',data.index(b'>',end))+1
    arrays={}
    for a in header.findall('.//PointData/DataArray'):
        pos=start+int(a.attrib['offset']);nbytes=struct.unpack_from('<Q',data,pos)[0]
        code='d' if a.attrib['type']=='Float64' else 'f'
        width=struct.calcsize(code)
        values=struct.unpack_from('<'+str(nbytes//width)+code,data,pos+8)
        arrays[a.attrib['Name']]=(int(a.attrib['NumberOfComponents']),values)
    return arrays


def legacy_arrays(path,binary):
    data=path.read_bytes();m=re.search(rb'\nPOINT_DATA (\d+)\n',data);assert m
    points=int(m[1]);rest=data[m.end():];arrays={}
    while rest.strip():
        rest=rest.lstrip(b'\r\n');line,rest=rest.split(b'\n',1);parts=line.split()
        assert parts[0] in (b'VECTORS',b'SCALARS'),line
        components=3 if parts[0]==b'VECTORS' else 1
        if components==1:
            lookup,rest=rest.split(b'\n',1);assert lookup==b'LOOKUP_TABLE default'
        count=points*components
        if binary:
            code='f' if parts[2]==b'float' else 'd';n=count*struct.calcsize(code)
            values=struct.unpack('>'+str(count)+code,rest[:n]);rest=rest[n:]
        else:
            tokens=rest.split(None,count);values=tuple(float(v) for v in tokens[:count])
            rest=tokens[count] if len(tokens)>count else b''
        arrays[parts[1].decode()]=(components,values)
    return arrays


def expected(dim,selection,ranks):
    complete=selection[0]=='Momentum_x' and selection[1:dim]==['Momentum_y','Momentum_z'][:dim-1]
    layout=[]
    if complete:
        layout.append(('Momentum',3,list(range(dim,2*dim))))
        for i,name in enumerate(selection[dim:],2*dim):layout.append((name,1,[i]))
    else:
        layout=[(name,1,[dim+i]) for i,name in enumerate(selection)]
    arrays={}
    for name,components,indices in layout:
        values=[]
        for rank in ranks:
            for point in range(8 if dim==3 else 4):
                values.extend(rank*1000+point*10+i for i in indices)
                if components==3 and dim==2:values.append(0)
        arrays[name]=(components,tuple(values))
    return arrays


def check(binary,ranks,mpi,launcher_name):
    with tempfile.TemporaryDirectory(prefix='output-writers-',dir=Path.cwd()) as tmp:
        work=Path(tmp)
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',TMPDIR=tmp,TMP=tmp,TEMP=tmp,
                 OMPI_MCA_orte_tmpdir_base=tmp,OMPI_MCA_btl_vader_backing_directory=tmp)
        launcher=[launcher_name,'-n',str(ranks)] if mpi else []
        def expected_error(result,messages):
            assert result.returncode!=0,result.stdout
            # An early aborting rank can prevent the MPI launcher forwarding rank zero's message.
            assert any(m in result.stdout for m in messages) or (ranks>1 and 'MPI_ABORT was invoked' in result.stdout),result.stdout
        def run(mode):
            folder=work/mode;folder.mkdir()
            result=subprocess.run(launcher+[str(binary),mode],cwd=folder,env=env,
                                  stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=60)
            (folder/'run.log').write_text(result.stdout)
            return folder,result
        folder,result=run('normal');assert result.returncode==0,result.stdout
        assert (folder/'accounting.bin').stat().st_size==5+3*ranks
        for name,reported in re.findall(r'WRITER_BYTES (\S+) ([\d.e+-]+)',result.stdout):
            path=folder/name;actual=path.stat().st_size
            if path.suffix=='.pvtu':actual+=sum(p.stat().st_size for p in path.with_suffix('').glob('*.vtu'))
            assert float(reported)==actual,(name,reported,actual)
        selections=[['Momentum_x'],['Momentum_y'],['Momentum_z'],['Momentum_x','Momentum_z'],
                    ['Momentum_x','Other_y','Momentum_z'],['Flux_x_extra'],
                    ['Momentum_x','Momentum_y','Momentum_z','Energy']]
        for dim in (2,3):
            for i,selection in enumerate(selections):
                if dim==2 and i==6:selection=[s for s in selection if s!='Momentum_z']
                prefix='d%ds%d'%(dim,i)
                full=expected(dim,selection,range(ranks))
                assert xml_arrays(folder/(prefix+'-xml.vtu'))==full,prefix
                assert legacy_arrays(folder/(prefix+'-ascii.vtk'),False)==full,prefix
                assert legacy_arrays(folder/(prefix+'-binary.vtk'),True)==full,prefix
                index=ET.parse(folder/(prefix+'-pieces.pvtu'))
                declarations={a.attrib['Name']:int(a.attrib['NumberOfComponents']) for a in index.findall('.//PPointData/PDataArray')}
                assert declarations=={name:v[0] for name,v in full.items()}
                for rank in range(ranks):
                    path=folder/(prefix+'-pieces')/(prefix+'-pieces_%d.vtu'%rank)
                    assert xml_arrays(path)==expected(dim,selection,[rank]),(prefix,rank)
        for rank in range(ranks):assert xml_arrays(folder/'empty'/('empty_%d.vtu'%rank))=={'Momentum_y':(1,())}
        for mode in ('limit-node','limit-product','limit-restart'):
            _,result=run(mode);expected_error(result,('32-bit',))
        if os.name=='posix':
            _,result=run('shared-error')
            messages=('MPI output file is incomplete','MPI output write failed','Closing MPI output failed','Output write failed','Closing output failed')
            expected_error(result,messages)
        if Path('/dev/full').exists() and os.name=='posix':
            folder,result=run('piece-error')
            expected_error(result,('Unable to write or close a PVTU piece',))
            assert not (folder/'failed-pieces.pvtu').exists()
    print('PASS: %d ranks, complete/partial 2D/3D vectors in four formats, empty pieces, I/O failure, size limits and byte/time accounting'%ranks)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--ranks',type=int,default=1);p.add_argument('--mpi',action='store_true');p.add_argument('--launcher',default='mpiexec')
    a=p.parse_args();check(a.binary.resolve(),a.ranks,a.mpi,a.launcher)
