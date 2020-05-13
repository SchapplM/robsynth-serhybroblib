% Calculate joint inertia matrix with Newton Euler for
% palh2m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:04
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1IC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1IC_inertiaJ_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1IC_inertiaJ_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1IC_inertiaJ_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1IC_inertiaJ_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1IC_inertiaJ_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:04:19
% EndTime: 2020-05-03 01:04:21
% DurationCPUTime: 1.58s
% Computational Cost: add. (1135->167), mult. (1634->195), div. (0->0), fcn. (782->8), ass. (0->93)
t1021 = (m(6) * pkin(6));
t908 = cos(qJ(5));
t1018 = -0.2e1 * t908;
t904 = sin(qJ(5));
t1019 = Ifges(6,4) * t1018 * t904 + Ifges(5,1) + Ifges(6,2) + (Ifges(6,1) - Ifges(6,2)) * t908 ^ 2 + ((2 * mrSges(6,3) + t1021) * pkin(6));
t985 = mrSges(6,2) * t904;
t993 = pkin(4) * m(6);
t1000 = -Ifges(5,2) - Ifges(6,3) + (mrSges(6,1) * t1018 + 0.2e1 * t985 - t993) * pkin(4) + t1019;
t909 = cos(qJ(4));
t897 = t909 ^ 2;
t1011 = t1000 * t897;
t861 = mrSges(6,1) * t904 + mrSges(6,2) * t908;
t1017 = mrSges(5,3) + t861;
t1020 = t1017 * pkin(3) - Ifges(4,5);
t843 = -(mrSges(6,2) * pkin(6) - Ifges(6,6)) * t904 + (mrSges(6,1) * pkin(6) - Ifges(6,5)) * t908;
t998 = mrSges(6,3) * pkin(4) + pkin(6) * t993 + Ifges(5,4) + t843;
t940 = 0.2e1 * t1000;
t933 = pkin(3) ^ 2;
t987 = t933 * m(6);
t877 = -mrSges(5,2) + mrSges(6,3) + t1021;
t905 = sin(qJ(4));
t865 = t877 * t905;
t863 = pkin(3) * t865;
t920 = m(5) * t933;
t1001 = Ifges(4,2) + t920 + 0.2e1 * t863;
t999 = -Ifges(4,1) + t1001;
t1016 = -t940 - 0.2e1 * t987 - 0.2e1 * t999;
t862 = mrSges(6,1) * t908 - t985;
t848 = mrSges(5,1) + t862 + t993;
t977 = t848 * t905;
t956 = pkin(3) * t977;
t1015 = -Ifges(4,4) + t956 + t998;
t914 = m(5) + m(6);
t1014 = pkin(3) * t914 + mrSges(4,1);
t847 = pkin(3) * t848;
t1013 = -0.2e1 * t847;
t989 = pkin(3) * t877;
t815 = t1000 * t905 + t989;
t1012 = t815 * t909;
t842 = t848 * t909;
t844 = t1014 + t865;
t820 = t842 + t844;
t906 = sin(qJ(3));
t910 = cos(qJ(3));
t1010 = -Ifges(4,6) * t910 + t1020 * t906;
t1008 = 0.2e1 * t820;
t1004 = t842 + t865;
t845 = 0.2e1 * t847;
t974 = t877 * t906;
t864 = pkin(2) * t974;
t1002 = -t864 + t1013 + t845;
t997 = Ifges(4,6) * t906 + t1020 * t910;
t996 = -0.8e1 * t998 * t905 - 0.4e1 * t847;
t994 = -0.2e1 * pkin(2);
t915 = m(4) + m(5);
t907 = sin(qJ(2));
t991 = pkin(1) * t907;
t990 = pkin(2) * t910;
t988 = t905 * pkin(3);
t981 = t914 * t933 + Ifges(4,3);
t980 = t998 * t897;
t839 = mrSges(4,2) + t977;
t979 = t839 * t906;
t978 = t844 * t906;
t976 = t848 * t906;
t972 = t906 * t907;
t971 = t906 * t909;
t948 = 0.4e1 * t998;
t970 = -t948 * t905 + t1013;
t934 = pkin(2) ^ 2;
t963 = t839 * t994;
t969 = t906 * t963 + t915 * t934;
t967 = pkin(3) * t862 + t843 * t905;
t965 = t861 * pkin(6);
t964 = -0.2e1 * t991;
t960 = -0.4e1 * t980;
t957 = pkin(2) * t976;
t952 = t905 * t957 + Ifges(4,3) + t863 + t920;
t825 = 0.2e1 * t998;
t854 = pkin(4) * t862 + Ifges(6,3);
t947 = t843 * t909 - t854 * t905;
t853 = t905 * t910 + t971;
t937 = t940 * t897 - t1000 - t999;
t935 = pkin(1) ^ 2;
t911 = cos(qJ(2));
t898 = t910 ^ 2;
t852 = t905 * t906 - t909 * t910;
t831 = -t957 + t989;
t819 = t852 * t907 - t853 * t911;
t818 = t852 * t911 + t853 * t907;
t807 = t1010 * t911 + t997 * t907;
t806 = (-Ifges(3,6) + t1010) * t911 + (-Ifges(3,5) + (mrSges(4,3) + t1017) * pkin(2) + t997) * t907;
t1 = [(t909 * t970 + t937 - t987) * t898 + (t906 * t960 + 0.2e1 * (-t815 * t906 + t877 * t991) * t909 + t839 * t964 + 0.2e1 * t1015 * t906) * t910 - t1011 + (t825 * t905 + t964 * t976 + t845) * t909 + (mrSges(3,2) + t978) * t964 + (t933 + t935) * m(6) + (m(3) + t915) * t935 + Ifges(3,1) + Ifges(2,3) + t1001 + (0.4e1 * t972 * t910 * t1011 + 0.2e1 * (t820 * t910 + t974 * t909 - t979 + mrSges(3,1) + (m(6) + t915) * pkin(2)) * pkin(1) + (0.4e1 * (t825 * t897 + t1012 - t1015) * t898 + ((0.2e1 * pkin(2) * t877 + t996 * t906) * t909 + t1016 * t906 + t963) * t910 + t960 + (-t905 * t940 - 0.2e1 * t957 - 0.2e1 * t989) * t909 + t978 * t994 + (2 * Ifges(3,4)) - 0.2e1 * Ifges(4,4) + t825 + 0.2e1 * t956) * t907 + (t937 + (0.2e1 * t864 + t970) * t909 + (-t933 + t934) * m(6) + (pkin(2) * t1008 + (0.4e1 * Ifges(4,4) - t948 - 0.4e1 * t956 + 0.8e1 * t980 + 0.4e1 * t1012) * t906) * t910 + (-t996 * t909 - 0.4e1 * t1011 - t1016) * t898 - Ifges(3,1) + Ifges(3,2) + t969) * t911) * t911 + t1019, t806, t807, ((t854 * t909 + t967) * t910 + pkin(2) * t862 + t947 * t906) * t911 - t967 * t972 + pkin(1) * t862 + (-t854 * t971 + t910 * t947) * t907; t806, -t831 * t905 + Ifges(3,3) + (t934 + t933) * m(6) + (t1008 - 0.2e1 * t1004) * t990 + t952 + t969, -t863 + t987 + (t864 + t1002) * t909 + (-t979 + (t820 - t1004) * t910) * pkin(2) + t952, (-pkin(2) * t853 - t988) * t861; t807, (-t831 + t989) * t905 + t1002 * t909 + (t1014 * t910 + (t877 * t909 - t839) * t906) * pkin(2) + t981, t981, -t861 * t988; -Ifges(6,3) * t818 + (Ifges(6,5) * t908 - Ifges(6,6) * t904) * t819 + t862 * (pkin(1) + t911 * pkin(2) - (-t910 * t911 + t972) * pkin(3) - t819 * pkin(6) - t818 * pkin(4)), -t861 * (pkin(2) * t971 + t905 * (pkin(3) + t990) + pkin(6)) + t965, -t861 * (pkin(6) + t988) + t965, Ifges(6,3);];
Mq = t1;
