import logging
import os
import shutil
import time
from abc import ABC, abstractmethod
from typing import List, Dict

import lcbuilder.constants
import requests
from lcbuilder.lcbuilder_class import LcBuilder
from astroquery.mast import Observations


class DvrUrlBuilder(ABC):
    """Base class to prepare the DVR url for a given mission id"""
    @abstractmethod
    def prepare(self, id: int, sectors: List[int], download_dir: str) -> str:
        """
        Prepares the DVR url for a given mission id
        :param int id: the id of the target in the given mission
        :param List[int] sectors: the sectors/quarters/campaigns of the DVR
        :param str download_dir: the directory to store the resulting files
        :return List[str]: The DVR paths
        """
        pass


class KeplerDvrUrlBuilder(DvrUrlBuilder):
    """Prepares the DVR url for a given Kepler mission id"""
    def prepare(self, id: int, sectors: List[int], download_dir: str) -> List[str]:
        """
        Prepares the DVR url for a given Kepler mission id
        :param id: the KIC id
        :param List[int] sectors: the quarters of the DVR
        :param str download_dir: the directory to store the resulting files
        :return List[str]: The DVR paths (only one for kepler targets)
        """
        ten_position_id = "{:10d}".format(id)
        three_first_numbers = ten_position_id[0:2]
        six_first_numbers = ten_position_id[0:5]
        filename = "dv/kplr{:s}-20160209194854_dvr.pdf".format(ten_position_id)
        url = "https://exoplanetarchive.ipac.caltech.edu/data/KeplerData/" \
              "{:s}/{:s}/{:s}/dv/kplr{:s}-20160209194854_dvr.pdf"\
            .format(three_first_numbers, six_first_numbers, ten_position_id, ten_position_id)
        tries: int = 1
        while tries < 4:
            try:
                response = requests.get(url)
                if response.status_code == 404:
                    logging.error("There is no official DVR for %s", id)
                    return []
                if response.status_code != 200:
                    raise ValueError("Unexpected status code.")
                with open(download_dir + '/' + filename, "wb") as f:
                    f.write(response.content)
                break
            except:
                logging.exception("Could not write the requested DVR for %s, retry no. %s",
                                  id, tries)
                time.sleep((tries + 1) ** 3)
            finally:
                tries = tries + 1
        return [download_dir + '/' + filename]


class TessDvrUrlBuilder(DvrUrlBuilder):
    """Prepares the DVR urls for a given Kepler mission id"""
    def prepare(self, id: int, sectors: List[int], download_dir: str) -> List[str]:
        """
        Prepares the DVR url for a given Kepler mission id
        :param int id: the KIC id
        :param List[int] sectors: the sectors of the DVR
        :param str download_dir: the directory to store the resulting files
        :return List[str]: The DVR paths
        """
        obsTable = Observations.query_criteria(provenance_name=[lcbuilder.constants.SPOC_AUTHOR,
                                                                lcbuilder.constants.TESS_SPOC_AUTHOR],
                                               target_name=[str(id)],
                                               sequence_number=sectors)
        data_products = Observations.get_product_list(obsTable)
        data_products = data_products[data_products['productSubGroupDescription'] == 'DVR']
        print(data_products['productFilename'])
        download_manifest = Observations.download_products(data_products, download_dir=download_dir)
        inner_download_dir = download_dir + '/mastDownload/TESS/'
        result_files = []
        for dir in sorted(os.listdir(inner_download_dir)):
            report_inner_dir = inner_download_dir + '/' + dir
            for file in sorted(os.listdir(report_inner_dir + '/')):
                if file.endswith('.pdf'):
                    result_path = download_dir + '/' + file
                    shutil.move(report_inner_dir + '/' + file, result_path)
                    result_files.append(result_path)
        shutil.rmtree(download_dir + '/mastDownload')
        return result_files


class DvrPreparer:
    """Class offering official data validation report download"""
    URL_BUILDERS: Dict[str, DvrUrlBuilder] = {lcbuilder.constants.MISSION_KEPLER: KeplerDvrUrlBuilder(),
                                              lcbuilder.constants.MISSION_TESS: TessDvrUrlBuilder()}
    """Dictionary containing the files preparer for each supported mission"""

    def retrieve(self, target_id: str, sectors: List[int], destination: str) -> List[str]:
        """
        Prepares the proper url for the given target_id data validation report and downloads it to the destination file
        :param str target_id: the target id
        :param List[int] sectors: the sectors of the DVR
        :param str destination: the destination file
        :return List[str]: list of downloaded files
        """
        try:
            mission, mission_prefix, id = LcBuilder().parse_object_info(target_id)
            logging.warning("Preparing official DVR for mission %s and target %s", mission, id)
            url_builder = self.URL_BUILDERS[mission]
            if url_builder is None:
                logging.warning("There is no official DVR for mission %s", mission)
            return self.URL_BUILDERS[mission].prepare(id, sectors, destination)
        except:
            logging.exception("Some problem when generating official DVR")
            return []
