% Calculate kinetic energy for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = picker2Dm1OL_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_energykin_fixb_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_energykin_fixb_slag_vp2: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_energykin_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1OL_energykin_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1OL_energykin_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:55
% EndTime: 2020-05-11 05:44:56
% DurationCPUTime: 0.15s
% Computational Cost: add. (187->58), mult. (326->106), div. (0->0), fcn. (136->14), ass. (0->39)
t138 = pkin(1) * qJD(1);
t120 = qJD(1) + qJD(2);
t128 = sin(qJ(2));
t137 = t128 * t138;
t118 = qJD(3) + t120;
t117 = qJD(4) + t120;
t134 = cos(qJ(2));
t115 = t134 * t138;
t111 = t120 * pkin(2) + t115;
t127 = sin(qJ(3));
t133 = cos(qJ(3));
t104 = -t133 * t111 + t127 * t137;
t110 = t120 * pkin(3) + t115;
t126 = sin(qJ(4));
t132 = cos(qJ(4));
t103 = t132 * t110 - t126 * t137;
t131 = cos(qJ(6));
t130 = cos(qJ(8));
t129 = cos(qJ(9));
t125 = sin(qJ(6));
t124 = sin(qJ(8));
t123 = sin(qJ(9));
t122 = cos(qJ(10));
t121 = sin(qJ(10));
t119 = qJD(1) + qJD(8);
t116 = qJD(6) + t120;
t114 = qJD(9) + t118;
t113 = qJD(10) + t117;
t108 = (-t125 * t134 - t128 * t131) * t138;
t107 = (t125 * t128 - t131 * t134) * t138;
t106 = -t127 * t111 - t133 * t137;
t105 = t126 * t110 + t132 * t137;
t102 = t118 * pkin(6) + t104;
t101 = t117 * pkin(4) + t103;
t100 = -t123 * t102 - t129 * t106;
t99 = -t129 * t102 + t123 * t106;
t98 = -t121 * t101 - t122 * t105;
t97 = -t122 * t101 + t121 * t105;
t1 = m(11) * (t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(10) * (t100 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t103 ^ 2 + t105 ^ 2) / 0.2e1 + m(4) * (t104 ^ 2 + t106 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2) / 0.2e1 + qJD(7) ^ 2 * Ifges(8,3) / 0.2e1 + qJD(5) ^ 2 * Ifges(6,3) / 0.2e1 + t119 ^ 2 * Ifges(9,3) / 0.2e1 + t120 ^ 2 * Ifges(3,3) / 0.2e1 + (t104 * mrSges(4,1) - t106 * mrSges(4,2) + Ifges(4,3) * t118 / 0.2e1) * t118 + (t103 * mrSges(5,1) - t105 * mrSges(5,2) + Ifges(5,3) * t117 / 0.2e1) * t117 + (t107 * mrSges(7,1) - t108 * mrSges(7,2) + Ifges(7,3) * t116 / 0.2e1) * t116 + (t99 * mrSges(10,1) - t100 * mrSges(10,2) + Ifges(10,3) * t114 / 0.2e1) * t114 + (t97 * mrSges(11,1) - t98 * mrSges(11,2) + Ifges(11,3) * t113 / 0.2e1) * t113 + ((mrSges(3,1) * t134 - mrSges(3,2) * t128) * t120 + (-mrSges(9,1) * t130 + mrSges(9,2) * t124) * t119) * t138 + (Ifges(2,3) / 0.2e1 + (m(9) * (t124 ^ 2 + t130 ^ 2) / 0.2e1 + m(3) * (t128 ^ 2 + t134 ^ 2) / 0.2e1) * pkin(1) ^ 2) * qJD(1) ^ 2;
T = t1;
