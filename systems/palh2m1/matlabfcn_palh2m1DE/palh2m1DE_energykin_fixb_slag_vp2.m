% Calculate kinetic energy for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m1DE_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1DE_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:51:57
% EndTime: 2020-05-02 23:51:57
% DurationCPUTime: 0.23s
% Computational Cost: add. (151->64), mult. (238->106), div. (0->0), fcn. (100->16), ass. (0->32)
t174 = qJD(1) / 0.2e1;
t156 = qJ(2) + qJ(3);
t160 = cos(qJ(2));
t167 = pkin(2) * qJD(2);
t150 = qJD(2) + qJD(3);
t168 = pkin(3) * t150;
t141 = cos(t156) * t168 + t160 * t167;
t173 = t141 ^ 2;
t172 = m(4) / 0.2e1;
t171 = m(5) / 0.2e1;
t170 = t174 - t150 / 0.2e1;
t169 = t174 - qJD(2) / 0.2e1;
t157 = sin(qJ(3));
t158 = sin(qJ(2));
t166 = t158 * t157;
t155 = -qJ(4) + qJ(2);
t154 = qJ(4) + qJ(2);
t165 = qJD(1) * (pkin(1) + pkin(4));
t159 = cos(qJ(3));
t152 = qJD(1) + qJD(2);
t151 = qJD(1) + qJD(4);
t149 = qJ(3) + t155;
t148 = qJ(3) + t154;
t146 = qJD(1) + t150;
t145 = t160 * pkin(2) + pkin(1);
t144 = (-t157 * t160 - t158 * t159) * qJD(1);
t143 = (-t159 * t160 + t166) * qJD(1);
t140 = -sin(t156) * t168 - t158 * t167;
t139 = (t159 * pkin(3) + pkin(2)) * t160 - pkin(3) * t166 + pkin(1);
t138 = cos(qJ(4)) * t165 + (t146 * cos(t149) / 0.2e1 + cos(t148) * t170) * pkin(3) + (t152 * cos(t155) / 0.2e1 + cos(t154) * t169) * pkin(2);
t137 = sin(qJ(4)) * t165 + (-t146 * sin(t149) / 0.2e1 + sin(t148) * t170) * pkin(3) + (-t152 * sin(t155) / 0.2e1 + sin(t154) * t169) * pkin(2);
t1 = (t140 ^ 2 + t173) * t171 + m(6) * (t137 ^ 2 + t138 ^ 2 + t173) / 0.2e1 + Ifges(4,3) * t150 ^ 2 / 0.2e1 + (Ifges(4,5) * t150 + Ifges(4,1) * t144 / 0.2e1) * t144 + (t138 * mrSges(6,1) - t137 * mrSges(6,2) + Ifges(6,3) * t151 / 0.2e1) * t151 + (Ifges(4,4) * t144 + Ifges(4,6) * t150 + Ifges(4,2) * t143 / 0.2e1) * t143 + ((t157 * (-t150 * mrSges(4,2) + t143 * mrSges(4,3)) + t159 * (t150 * mrSges(4,1) - t144 * mrSges(4,3))) * pkin(2) + ((t157 ^ 2 + t159 ^ 2) * pkin(2) ^ 2 * t172 + Ifges(3,3) / 0.2e1) * qJD(2)) * qJD(2) + (t145 * (-t143 * mrSges(4,1) + t144 * mrSges(4,2)) - t140 * mrSges(5,3) + qJD(2) * (-Ifges(3,5) * t158 - Ifges(3,6) * t160) + (Ifges(2,3) / 0.2e1 + t145 ^ 2 * t172 + m(3) * pkin(1) ^ 2 / 0.2e1 + Ifges(5,2) / 0.2e1 + (t139 * t171 + mrSges(5,1)) * t139 + (pkin(1) * mrSges(3,1) + Ifges(3,2) * t160 / 0.2e1) * t160 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t160 + Ifges(3,1) * t158 / 0.2e1) * t158) * qJD(1)) * qJD(1);
T = t1;
