"""查看首个GRIB消息前导，不扫描场或解码打包数据。"""
from .io import asset

def summarize(v,root):
    raw,digest=asset(root,v['file'])
    if len(raw)<16 or raw[:4]!=b'GRIB' or raw[7]!=2:raise ValueError('首消息前导无效')
    return {'file_bytes':len(raw),'first_declared_message_bytes':int.from_bytes(raw[8:16],'big'),'edition':raw[7],'sha256':digest}
