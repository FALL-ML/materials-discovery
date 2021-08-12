import os

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#系统用到的 mongdb数据库设置
MONGODB = {
    'MaterialsDB': {
        'DatabaseName': 'materialDB',
        'User': '',
        'Password': '',
        'Host': '***',
        'Port': ****,
    },
        'Firworks': {
        'DatabaseName': 'fireworks',
        'User': '',
        'Password': '',
        'Host': '***',
        'Port': ****,
    }
}

#计算服务器配置信息
COMPUTER_SERVER = {
    "ziqiang4000":
    {
        "user": '***',
        "password": "",
        "port": ***,
        "hostaddr": '***',
        "hostname": "***",
        "exec_dir": '****',
    }
}
