% Calculate kinetic energy for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnOL_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:46
% EndTime: 2020-04-12 19:40:46
% DurationCPUTime: 0.13s
% Computational Cost: add. (57->33), mult. (155->75), div. (0->0), fcn. (72->6), ass. (0->11)
t51 = qJD(1) ^ 2;
t49 = cos(qJ(2));
t48 = cos(qJ(3));
t47 = cos(qJ(4));
t46 = sin(qJ(2));
t45 = sin(qJ(3));
t44 = sin(qJ(4));
t43 = qJD(2) + qJD(3);
t42 = (-t45 * t49 - t46 * t48) * qJD(1);
t41 = (t45 * t46 - t48 * t49) * qJD(1);
t1 = Ifges(4,3) * t43 ^ 2 / 0.2e1 + Ifges(5,3) * qJD(4) ^ 2 / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(5) * pkin(1) ^ 2 / 0.2e1) * t51 + (Ifges(4,5) * t43 + Ifges(4,1) * t42 / 0.2e1) * t42 + m(4) * (t49 ^ 2 * t51 + (t45 ^ 2 + t48 ^ 2) * qJD(2) ^ 2) * pkin(2) ^ 2 / 0.2e1 + (Ifges(4,4) * t42 + Ifges(4,6) * t43 + Ifges(4,2) * t41 / 0.2e1) * t41 + (Ifges(3,3) * qJD(2) / 0.2e1 + (-t48 * (t43 * mrSges(4,1) - t42 * mrSges(4,3)) - t45 * (-t43 * mrSges(4,2) + t41 * mrSges(4,3))) * pkin(2)) * qJD(2) + (t46 * Ifges(3,5) * qJD(2) + (Ifges(3,1) * t46 ^ 2 / 0.2e1 + (pkin(1) * mrSges(5,1) + Ifges(5,2) * t47 / 0.2e1) * t47 + (-pkin(1) * mrSges(5,2) + Ifges(5,4) * t47 + Ifges(5,1) * t44 / 0.2e1) * t44) * qJD(1) + qJD(4) * (Ifges(5,5) * t44 + Ifges(5,6) * t47) + (Ifges(3,6) * qJD(2) - pkin(2) * (-t41 * mrSges(4,1) + t42 * mrSges(4,2)) + (Ifges(3,4) * t46 + Ifges(3,2) * t49 / 0.2e1) * qJD(1)) * t49) * qJD(1);
T = t1;
