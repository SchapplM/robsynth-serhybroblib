% Calculate kinetic energy for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1DE1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_energykin_fixb_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_energykin_fixb_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_energykin_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE1_energykin_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1DE1_energykin_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1DE1_energykin_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:55:34
% EndTime: 2020-04-24 19:55:35
% DurationCPUTime: 0.06s
% Computational Cost: add. (220->24), mult. (321->39), div. (14->5), fcn. (74->4), ass. (0->20)
t135 = pkin(3) ^ 2;
t130 = pkin(1) * cos(qJ(1));
t128 = pkin(1) ^ 2 - 0.2e1 * pkin(2) * t130;
t132 = -pkin(3) + pkin(4);
t133 = -pkin(3) - pkin(4);
t134 = ((pkin(2) - t133) * (pkin(2) + t133) + t128) * ((pkin(2) - t132) * (pkin(2) + t132) + t128);
t131 = pkin(1) * sin(qJ(1));
t129 = pkin(2) * qJD(1);
t127 = -pkin(4) ^ 2 + t135;
t112 = pkin(2) ^ 2 + t128;
t126 = t129 * t131;
t123 = sqrt(-t134);
t113 = -pkin(2) + t130;
t110 = t112 + t127;
t103 = t113 * t123 * t129;
t102 = t103 + (t112 - t127) * t126;
t101 = -t110 * t126 + t103;
t100 = (t110 * t113 + t123 * t131) * t129;
t99 = t101 ^ 2;
t1 = qJD(1) ^ 2 * Ifges(2,3) / 0.2e1 + (m(3) * (t99 / 0.4e1 + t100 ^ 2 / 0.4e1) / t135 / 0.2e1 + (t101 * mrSges(3,2) / 0.2e1 - t100 * mrSges(3,1) / 0.2e1) * t102 / t123 / pkin(3) - (t102 ^ 2 * Ifges(3,3) + t99 * Ifges(4,3)) / t134 / 0.2e1) / t112 ^ 2;
T = t1;
