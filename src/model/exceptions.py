# Copyright 2018 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


class NoIDMapping(Exception):
    """Thrown when a search for a given metabolite identifier does not yield any result."""
    def __init__(self, compound_id):
        self.value = compound_id

    def __str__(self):
        return f"No metabolite associated with {self.value}"


class PartNotFound(Exception):
    """Thrown by the ICE client when a requested genotype part is not found"""
    pass
