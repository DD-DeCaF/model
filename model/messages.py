from venom.fields import String, Repeat
from venom.message import Message
from venom.rpc import Stub
from venom.rpc.stub import RPC


class GeneRequest(Message):
    gene = String()


class ReactionsResponse(Message):  # use Map() field when implemented
    reactions_ids = Repeat(String())
    equations = Repeat(String())


class GeneToReactionsRemote(Stub):
    reactions = RPC(GeneRequest, ReactionsResponse)
